# Methods (working notes)

These notes record how the detector works and why each piece is there. Numbers quoted are
from the development case (Denver, 25-26 January 2025) unless stated otherwise.

## Data

| Source | What | Notes |
|---|---|---|
| adsb.lol history (ODbL) | Position of every tracked aircraft, one per 10 s slot | Heatmap section of the daily archive; pressure altitude in 25 ft steps (GNSS altitude if baro is missing); kept below 30 000 ft over the US and Alaska |
| NEXRAD Level III N0B (AWS) | 0.5 deg super-resolution reflectivity, 720 rays x 0.25 km gates | Gate-for-gate identical to the Level II lowest sweep (checked) |
| HRRR analyses (Utah Zarr archive) | Hourly temperature, humidity, cloud water and wind on 25 hPa levels | Grid verified against the published corner point to < 1 m |
| tar1090-db | Aircraft address -> registration, ICAO type, engine class | ICAO type description gives jet / turboprop / piston / helicopter |
| IEM | METARs | Ceilings, weather, surface temperature |

The Denver radiosonde archive ends on 8 July 2022 (IGRA2 and IEM both), so vertical
profiles come from HRRR.

## Geometry

Radar gates are placed with the 4/3 effective earth radius model. Aircraft pressure
altitude converts directly to pressure with the ISA formula, which lets track points be
matched to HRRR pressure levels without any altimeter correction.

## The physical signature

Ice produced by an aircraft in supercooled cloud grows, falls and drifts with the wind. On
the lowest radar sweep (here 300-700 m above ground near the airport) it appears 10-40
minutes later as a narrow line of snow echo (1-3 km wide, 10-30 km long, 5-25 dBZ) that
follows the shape of the flight path, displaced downwind. In the Denver case the drift that
best aligned tracks with echo (0.5-1.0 m/s east, 3.0-4.5 m/s north) agreed with the HRRR
mean wind between the radar beam and the seeding level (0-0.3 east, 2.9-3.3 north).

## Aircraft echoes

Aircraft themselves are radar targets. Every tracked aircraft (any altitude) is masked
where the beam actually saw it: the sweep's start time, starting azimuth and duration (from
the Level III header, checked against Level II ray times: 71 s vs 71.8 s) give the time
the beam pointed at each azimuth; each aircraft's track is interpolated to that moment and a
disk of 1.5 km plus 5 s of flight is removed. This removes ~4% of the radar area at Denver.

## Test 1: parcel-following difference-in-differences (event.py)

Every track point's air parcel is followed with the drift from 30 min before to 60 min after
passage. For points grouped by altitude band (2000 ft, sliding by 1000 ft) and passage hour,
S = (change of reflectivity on the track) - (change on lines 6 km either side), between
"before" (-30 to -5 min) and "after" (15 to 60 min). Placebos are the same point sets
rotated about the radar by 20-340 degrees (161 rotations). Cells are converted to z-scores
against the placebos and the day is scored with a cluster-mass statistic (cells with z > 2
connected in band x hour), compared with the largest cluster of each placebo
(leave-one-out). Denver: p = 0.025.

Weakness found: near a busy airport the "side" lines are not empty - they fall on other
flights' paths and streaks - and the average over all points is diluted by the many
flights that do not cross cloud. Worse, it is not calibrated on winter days: on eight
January 2025 Denver days on which the streak test found no evidence of seeding (1-8
January, day p-values 0.08-0.99), its primary statistic reached p <= 0.05 on three (0.025,
0.012, 0.006), where 0.4 such days would be expected by chance. Real traffic is concentrated on a few corridors, and the change in echo along
those corridors is not exchangeable with that along rotated copies on days when
precipitation is not uniform. It is therefore not used in the census
(pipeline.run_day(event_test=False)); the code remains for reference.

## Test 2: streak objects matched to single flights (streaks.py, match.py)

1. Streak objects are found in each scan without reference to aircraft: Gaussian smoothing
   (0.75 km), white top-hat with a 4 km disk (removes broad precipitation), threshold
   (top-hat >= 6 dB and >= 3 dBZ), cleanup, then shape by skeleton: main path >= 8 km,
   mean width <= 3.5 km, main path >= 4 x mean width, branchiness <= 2.5, and peak
   reflectivity <= 35 dBZ (aircraft snow is light; brighter lines are convective).
   Skeleton measures accept the U-shaped streaks left by approach patterns.
2. Fixed echoes are removed: objects lying mostly on pixels that are inside streak objects
   in >= 40% of the day's scans (ground clutter on terrain edges).
3. Each object is compared with every flight that passed 10-50 min earlier. Tracks are
   interpolated to 250 m spacing and drifted with the HRRR drift plus a perturbation of
   up to +-1 m/s (0.5 m/s steps). Coverage (share of the object's centre line within 1 km
   of the drifted track) and precision (share of the drifted track near the object that is
   on it) are combined as F1. Points slower than 80 kt are not used: helicopters and slow
   aircraft loop over small areas, so their drifted "lines" are blobs that fit any echo
   (one such helicopter produced the only false-alarm day during development).
   Only physically possible seeding points are candidates (match.py, rotation-invariant so
   real and placebo flights are filtered alike):
   - fall speed: ice made dh above the lowest beam's centre (geometric height from the HRRR
     geopotential at the aircraft's pressure, minus the 4/3-earth beam height; beam
     half-width w = range x tan 0.475 deg) must reach the beam at the point's age a with a
     fall speed of 0.2-2.0 m/s (snow crystals to lightly rimed aggregates, Locatelli &
     Hobbs 1974): 0.2 a <= dh + w and dh - w <= 2.0 a;
   - temperature: HRRR temperature at the point <= -4 C (supercooled cloud needed).
   Each match reports its implied fall speed, median dh / a over its seeding points.
4. Newness: the object, moved back with the drift, must have been echo-free (>= 70% of
   samples < 0 dBZ) in the scans 3-30 min before the flight passed. In a stream of aircraft
   the match is therefore the first aircraft through new cloud.
5. Chance: the whole procedure is repeated with all flights rotated about the radar
   (24 rotations to screen a day, 200 to confirm a day that screens at p <= 0.05). A day is
   scored by E = sum over matched flights of -log10 of each flight's best streak p-value (a
   flight matched to the same echo in many scans counts once). Its null distribution treats
   each rotation in turn as the real set. The rank p-value cannot go below 1/(R + 1); a
   gamma distribution fitted to the null scores gives a smooth tail estimate for ranking
   days across the census (its calibration is checked on null days).
6. False discovery rate within a day: for a p-value threshold, the expected number of
   falsely matched flights is the mean number each rotation produces when treated as real.
7. Time-reversed diagnostic: the same matching with aircraft that pass 10-50 min after each
   scan (no newness test). They cannot have made the echo, so an excess over rotations would
   reveal streak-like echoes aligned with traffic corridors for other reasons. It is not a
   pure null: on a day with real seeding, the aircraft following a seeder along the same
   path match its streak in reverse.
8. Orientation check (scripts/shift_null.py): rotations change a track's orientation, so if
   streak-like echo shared the preferred orientation of real paths (bands and runways both
   follow the wind) rotations would favour the real flights. Translating every flight by
   8-16 km in 8 directions keeps orientation exactly; both placebo families are compared on
   the same tracks (within 49 km of the radar, so shifted copies stay in the domain).

Census selection (scripts/census_summary.py): every radar-day keeps its 24-rotation
screening score. The gamma-tail p-value is calibrated empirically: each rotation, scored
against the other rotations of its day, gives a null p-value; pooled over all tested days
these give the calibrated p-value of each real day, (1 + count of pooled null p-values at
or below it) / (1 + their number). Event days are selected by Benjamini-Hochberg at FDR 5%
on the calibrated p-values; a day that screens at rank p <= 0.05 or gamma p <= 0.1 is
re-run with 200 rotations for its per-flight catalogue and false discovery rates.

Denver, 25-26 January 2025 (analysis day from 12 UTC, final method): 677 streak objects,
63 matched to a flight; E = 78 against a null median of 29 and 95th percentile of 48; no
rotation of 200 reached it (rank p < 0.005, gamma tail p = 2e-4). At p <= 0.005, 14 flights
(16 streaks), all jets, with 0.7 expected false (FDR 5%).

## The Denver case in physical terms

- Seeding flights: CRJ-200, E175, A319/A320, 737-800/MAX 8, 757, 787-9, Citation Latitude;
  mostly arrivals; seeding altitudes 6 800-13 200 ft (pressure altitude); HRRR temperature at
  the seeding points -12 to -18.5 C, mostly -16 to -18.5 C.
- Streaks appeared at the lowest tilt 13-46 min after passage, were 9-30 km long and ~2-3 km
  wide, with peak reflectivity 14-29 dBZ (roughly 0.5-2.7 mm/h liquid equivalent with
  Ze = 110 S^2, Wolfe & Snider 2012; uncertain by a factor of a few).
- Composite over the matched flights: the drifted track stands 7-10 dB above lines 5-8 km
  to either side from ~17 to ~52 min after passage and is back to background by ~60 min.
  (The 10-50 min matching window and the newness test shape the composite's onset, so only
  its width and decay are free of selection.)
- Dual polarisation (Level II, lowest tilt): streak snow has lower differential
  reflectivity than other echo of the same reflectivity (median ZDR 0.22 vs 0.44 dB at
  10-20 dBZ; correlation coefficient ~1.00 vs 0.99), i.e. rounder, aggregated particles,
  in line with Moisseev et al. (2019) at Helsinki. Gates are not independent; per-streak
  medians are the fair comparison.
- Efficiency: of 48 jets that first crossed the -17 to -20 C layer (air echo-free before
  they arrived), 62% left a streak along their drifted track against 10% for the same
  flights rotated about the radar. Most flights in warmer bins flew below the cloud base, so
  this is not yet a temperature dependence.
- Surface reports (Buckley SFB, DEN, Centennial) show light snow on and off all evening, so
  the aircraft contribution cannot be separated at the ground from hourly METARs.

## Calibration

- Recall from synthetic streaks injected into real radar on three quiet Denver days (1-3
  January 2025), recovery = matched to the injected flight at p <= 0.05 against 101
  rotations. The injection is physical (inject.py, scripts/injection_recovery.py): a
  supercooled layer 300-800 m deep is placed at a height a real flight crossed; its points
  in the layer seed; each point's ice reaches the beam after (height above the beam) /
  (fall speed, 0.4-1.2 m/s, narrower than the detector's allowed range); nothing is
  visible before 10 min; the streak is 1.5-3.0 km wide above 3 dBZ when it appears (the
  widths of the detected Denver streaks) and widens at 0.2 m/s; HRRR wind error 0.7 m/s per
  component. Recall (138 trials, 95% intervals): 10/33 at 5 dBZ peak (30%, 16-49%), 15/30 at
  10 dBZ (50%, 31-69%), 19/34 at 15 dBZ (56%, 38-73%), 27/41 at 20 dBZ (66%, 49-80%). Most
  trials were flights climbing or descending through the layer (88 of 105 at >= 10 dBZ
  crossed more than 1000 ft of it), which leave shorter, sheared streaks: 55% of those
  were recovered, against 77% of flights within 300-1000 ft and 3 of 4 near-level ones.
- Two earlier injection models were wrong in instructive ways. Injecting every point of a
  6 000-14 000 ft band on a fixed timetable (onset 15, peak 35 min) gave 85-95% recall, but
  ignored that ice from different heights reaches the beam at different times. Adding that
  timing with a fixed Gaussian width made brighter streaks harder to find (42% found at
  20 dBZ against 73% at 10 dBZ): a Gaussian in linear units is wider above a fixed
  threshold the brighter it is, so bright injections exceeded the 3.5 km width limit, unlike
  real streaks.
- False alarms: the same pipeline on July days restricted to tracks below 13 000 ft, where
  the air is above freezing and aircraft cannot seed (null mode; only points the HRRR puts
  above 0 C are used).
- The rotation null fails on at least one warm day. On 1 July 2025 at Denver, with no
  seeding possible, the real flights matched streak-like echo better than all 200 rotated
  copies before the physical requirements were added (p = 0.005; 2 035 objects, 7 flights at
  p <= 0.005 against 1.6 expected), and better than all but one after them (p = 0.010;
  1 452 objects with peak <= 35 dBZ, 7 flights at p <= 0.005 against 0.9 expected; gamma-tail
  p = 0.006). The matched objects were edges of convective echo 20-35 km south-west of the
  radar, along the southern arrival corridor, matched by low flights (implied fall speeds
  0.12-0.6 m/s). Rotation breaks every link between where aircraft fly and where echo
  forms, and not all such links are seeding: aircraft are routed around storms, and
  terrain-forced convergence lines (the Denver Convergence-Vorticity Zone, for one) form near
  the airport. One warm day at p = 0.01 could be chance; the warm-season null at all census
  radars (below) measures how often it happens.
- A corridor-preserving null does not work, because airline schedules repeat. The day-swap
  test (scripts/dayswap_null.py) replaces the real flights with those of each of the 28
  other days within two weeks, at the same times of day, through this day's radar and wind.
  On 1 July it removes the false alarm (rank p = 0.33). But it also removes Denver (E = 31.8
  against a median of 22.3 and 95th percentile of 31.9 for the other days; rank p = 0.10):
  on 7 of the 12 confirmed streaks, flights of 6-14 of the 28 other days fit at least as
  well as the real seeder (F1 up to 0.94), because the same flight numbers fly the same
  procedures at the same minutes every day. A null that keeps the corridors and the timing
  keeps the clones of the seeders too, so it can say that a streak lies on a corridor but
  not that one aircraft made it. The census therefore uses rotations, with the physical
  requirements, and measures their false-alarm rate on warm days instead.

## Corrections after an independent code review

A read-only review of the statistics and pipeline code found these, all fixed:

- **Near-copy rotations made the null too easy to beat.** Rotations 1.6 deg apart move a
  track at 50 km by 1.4 km, less than a match can absorb (up to ~5 km with the drift
  correction), so neighbouring rotations match the same streaks. The real flights have no
  rotation within 20 deg of them; but a rotation scored "as if real" was compared with its
  near-copies, which deflated its evidence and so the null distribution of the day score.
  In a simulation with this structure the 200-rotation day p-value was <= 0.01 on 6.3% of
  null days. Now each set is compared only with sets at least 20 deg away, as the real one
  is (match.comparable_sets); that brings it to 4.7%, and the gamma tail to 3.7%. The rest
  is not a bug but a limit: there are only some 15-25 independent rotations, so no rotation
  p-value below about 1/20 is exact. Census selection therefore uses an empirical null
  (below), not the rotation p-value's tail.
- F1 ties: placebo scores were stored in float32 and the real score compared in float64, so
  ties broke in the real set's favour. All comparisons are now in float32.
- Injection recovery counted a match to any flight as a recovery; it now requires the
  injected flight.
- Smaller: unknown ground speed (-1) was interpolated into fake slow points that were
  dropped; the newness test counted scans that had not seen the air; points without a drift
  could crash matched_points; a missing gamma tail became p = 1 instead of the rank p-value.

## Non-meteorological echo

Two of the warm-null false alarms came from echo that is not precipitation. At Buffalo,
four of the nine significant matches on 1 July 2025 were at one spot 33 km south of the
radar, which held a streak-shaped object in 23 scans through the day; on 3 November 2024 at
Indianapolis the matched "streaks" were patches of speckled, radially streaked clear-air
echo. The correlation coefficient separates these: snow and rain are above ~0.97, ground
clutter, wind turbines, birds, insects and anomalous propagation well below 0.9. Each
reflectivity scan (Level III N0B) is now paired with the same scan's correlation
coefficient (N0C), and echo below 0.90 counts as no echo, before streak detection and the
newness test. The Denver streaks measured 0.992-0.995.

On the three warm null days that had come out significant, the mask removed most streak
objects and every false alarm: Denver 1 452 -> 60 objects, p = 0.010 -> 0.72; Buffalo 455 ->
73, p = 0.005 -> 0.76; Cleveland 729 -> 249, p = 0.010 -> 0.60. What had looked like edges
of convective echo near Denver was mostly clear-air echo: summer boundary layers are full of
insects, which radars see along convergence lines as thin bright lines with a low
correlation coefficient.

## Lessons from development

- A global average over all flights and all of the radar domain dilutes the signal to the
  point of uselessness; the signal is local because supercooled cloud is patchy.
- "Control" lines beside a track are not clean near busy airports.
- Time-shift placebos are contaminated by the real effect when the same paths are flown
  all day; rotation placebos are not.
- With many flights, many drift perturbations and long age windows, some random flight
  always fits; restricting ages to 10-50 min and requiring precision as well as coverage
  are what make the single-flight match specific.
- Heatmap positions are 2-6 km apart at airliner speeds; tracks must be interpolated before
  any geometric comparison at 1 km tolerance.
- Counting every streak object in every scan lets one persistent coincidence count many
  times; scoring flights, not objects, fixed the helicopter false alarm.
- numpy.isclose on epoch seconds is not an equality test (rtol 1e-5 x 1.7e9 s = 4.7 h);
  scan matching uses an absolute tolerance.
- Many airliners broadcast an outside air temperature, but it is derived from Mach and true
  airspeed and is unreliable below Mach ~0.5, so it cannot check HRRR at approach speeds.
- The Denver radiosonde archive ends in July 2022, so HRRR is the only profile source there.
