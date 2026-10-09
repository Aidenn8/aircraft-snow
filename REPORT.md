# Aircraft-made snow: an automated census from open radar and flight data

## Summary

Aircraft climbing or descending through clouds of supercooled droplets can freeze some of
them; the ice grows, falls and drifts with the wind, and weather radar sees a narrow band
of snow that traces the flight path, displaced downwind. This had been shown in individual
cases. This project builds an automated detector for these streaks in NEXRAD radar data,
matches each streak to the single aircraft that made it using open ADS-B tracks, and
measures how often chance alone would produce such a match.

On 25 January 2025 at Denver International Airport, the real flights' streak matches beat
those of every one of 200 copies of the day's flights rotated to other bearings
(gamma-tail p = 0.0013), while aircraft passing after each scan, which cannot have made its
echoes, match no better than chance. Sixteen flights matched streaks along their own
wind-carried tracks at p <= 0.01, against 3.6 expected by chance; their ice implied fall
speeds of 0.3-1.3 m/s, as for snow. About one jet in eight that crossed the supercooled
cloud layer into clear air left a radar-visible streak, nine times the chance rate.

Run without supervision on every winter day at eight radars near large airports (2,888
tested radar-days, November 2023 to March 2026 plus March 2023), the test found streaks
beyond chance on 58 days at p <= 0.01 where 29 would be false alarms, a rate measured on 744
warm days on which no aircraft can make ice (7 false alarms, 7.4 expected). The excess
comes on days with supercooled cloud over the radar and none on days without it, a
dependence that was predicted before the last winter was analysed and held there. A
catalogue of 163 flights on days judged convincing by eye (seven radars, 155 of them jets,
from regional jets to 747 freighters) shows seeding at a median of 7 300 ft and -13 C, much
warmer than usually quoted for jets, and arriving jets leaving streaks at about five times
the odds of departing ones in the same cloud.

## 1. Background

See `docs/background.md`. In short: aircraft-produced ice particles (APIPs) were shown in
field tests in the 1980s-2000s; Heymsfield et al. (2010, 2011) documented hole-punch and
canal clouds from commercial aircraft and estimated that conditions allowing it occur up
to ~6% of the time at major airports; Moisseev et al. (2019) found aircraft-intensified snow
by inspecting 11 years of Helsinki radar; Iwasaki & Fujii (2025) matched 59 hole-punch clouds
in satellite imagery to Flightradar24 tracks. No study has counted radar-visible
aircraft-made snow across airports and winters, or measured which aircraft make it and at
what temperatures.

## 2. Data

| Source | Use |
|---|---|
| adsb.lol open ADS-B history (ODbL), March 2023 onward | aircraft positions every 10-20 s, pressure altitude |
| NEXRAD Level III N0B and N0C (AWS) | lowest-tilt reflectivity (0.5 deg x 250 m) and correlation coefficient |
| NEXRAD Level II (AWS) | dual-polarisation moments for the case study |
| HRRR analyses (Utah Zarr archive) | hourly temperature, humidity, wind profiles |
| tar1090-db | aircraft type and engine class |
| Iowa Environmental Mesonet | METARs |

## 3. Method

See `docs/methods.md` for the full description and the lessons from building it. In brief:

1. Echo that is not precipitation (insects, birds, ground clutter, wind turbines; a
   correlation coefficient below 0.90) is removed, and streak-shaped echo objects (long,
   thin, at most 35 dBZ) are found in each radar scan by morphology alone.
2. Each is compared with every flight that passed 10-50 minutes earlier, its track carried
   by the HRRR wind between the radar beam and the flight's level (with up to 1 m/s of
   allowed wind error). The match must cover the streak and lie along it, and the air must
   have been echo-free before the aircraft arrived.
3. Only physically possible seeding counts: the aircraft must have been in air at or below
   -4 C, and high enough above the radar beam that ice falling at 0.2-2 m/s (the range for
   snow) would reach the beam at the observed age.
4. Chance is measured by repeating everything with all flights rotated about the radar,
   and checked with flights translated instead (which keeps their orientation to the wind).
5. A radar-day is scored by the combined evidence of its matched flights, against the same
   score for each rotated flight set; across the census, days are selected at a 5% false
   discovery rate after calibrating the score's tail on the rotations themselves.

Aircraft echoes are removed from the radar where the beam saw each aircraft, helicopters and
slow aircraft are excluded from matching, and fixed clutter is removed by persistence.

## 4. The Denver case, 25 January 2025

![Six matched aircraft and their streaks](figures/denver_case.png)

*Lowest-tilt reflectivity of the Denver radar (KFTG) at six times, with the traffic of the
previous 45 minutes (grey) and, in orange, the air each matched aircraft flew through,
carried by the wind to the scan time.*

- **Conditions.** An easterly upslope flow (6-11 kt at the airport, -8 C at the surface)
  under broken to overcast cloud at 3 500-5 000 ft above the airport (mostly 4 600-4 800 ft,
  about 3.1 km above sea level), light snow on and off. HRRR has a moist layer (83-88%
  relative humidity) up to about 3.3 km at -18 to -19 C, capped by an inversion that warms
  to -12 C near 3.9 km; it has no cloud water in it, a known weakness for thin supercooled
  stratus.
- **The day as a whole.** Of 685 streak-shaped echoes, 65 were matched to a flight. The
  day's combined evidence (E = 78.7) exceeded that of every one of 200 copies of the day's
  flights rotated to other bearings (median 34, 95th percentile 53; gamma-tail p = 0.0013).
  Aircraft passing 10-50 minutes *after* each scan, which cannot have made its echoes,
  match no better than chance (p = 0.12).
- **Which aircraft.** At p <= 0.01, 16 flights matched streaks along their own
  wind-carried tracks, against 3.6 expected by chance: about twelve genuine seeders, with
  roughly one named match in four possibly chance. Fifteen were jets (CRJ-200 x3, E175 x2,
  737-800 x2, 737 MAX 8 x2, A321neo x2, 737-900, 787-9, A319, Citation Latitude) and one a
  Pilatus PC-12 turboprop; 11 were arriving, 4 departing. Nine of the 16 crossed the cloud
  at 9 700-10 900 ft (about 3.0 km, near the reported cloud base), and 11 at HRRR
  temperatures of -16 to -18.5 C.
  Two matches look physically odd and are the likeliest chance ones: an A319 seeding at
  20 000 ft and -26 C, and an A321neo only 500 m above the radar beam; the PC-12 match, at
  dawn in widespread snow, is also doubtful on inspection.
- **Streaks.** The matched streaks (20, from 13:34 to 05:16 UTC) were 9-28 km long and
  1.4-3.3 km wide (median 2.3 km), peaking at 13.5-29.5 dBZ, and were seen 15-49 minutes
  after the aircraft passed; most lie along the downwind, base and final legs of the arrival
  pattern south-west of the airport.
- **Physical consistency.** The seeding points were 0.5-3.7 km (median 0.9 km) above the
  radar beam, which over the observed delays implies fall speeds of 0.3-1.3 m/s (median
  0.55 m/s; 18 of 20 within 0.35-1.0 m/s), the speeds of snow crystals and aggregates.
  Allowing the wind to be wrong by up to 3 m/s instead of 1 m/s improves the fits by a
  median of 0.001 in F1.
- **Life cycle.** Following each matched flight's air through every scan, echo on the
  drifted track is no higher than beside it before the aircraft passes, rises within
  minutes, peaks at about 8.5 dB above lines 5-8 km either side 25-35 minutes later, and is
  back to background after about an hour. (Between 10 and 50 minutes the composite is
  shaped by how flights were selected; its flat start and its decay are not.)

  ![Composite life cycle](figures/denver_lifecycle.png)

- **Particles.** Compared scan by scan with other echo of the same strength, streak echo
  has a higher correlation coefficient (15 of 17 scans at 10-20 dBZ, 10 of 10 at 20-30 dBZ;
  medians 0.995 against 0.988 and 0.974), i.e. more uniform particles. Its differential
  reflectivity is lower only for weak echo (0-10 dBZ, 14 of 17 scans), which may say as
  much about the noisier background echo.
- **How many jets did it.** Of the jets that entered echo-free air in the seeding layer
  (8 650-11 800 ft, 21:20-04:40 UTC; first passages), 13.6% left a streak along their
  drifted track, against 1.5% for the same flights rotated about the radar: about one jet
  in eight, nine times the chance rate. The share was 15.7% at -17 to -20 C (53 of 338) and
  11.0% at -14 to -17 C (31 of 282); the difference is within the uncertainty and could come
  from height above the beam as much as temperature. Ten turboprops crossed; one left a
  streak.

  ![Efficiency by temperature](figures/denver_efficiency.png)

- **How much snow.** With Ze = 110 S^2 (Wolfe & Snider 2012), the streaks peaked at
  0.5-2.9 mm/h of snow (water equivalent), and each was depositing roughly 4 000-38 000
  tonnes of water per hour at its matched scan (median 15 000). Over half an hour that is
  comparable to all the liquid water a cloud this thin holds, so either the ice also grew
  from vapour or the relation overstates snowfall for young, small crystals; the amounts
  are uncertain by a factor of several.
- **Ground.** Stations at Denver International, Buckley SFB and Centennial reported light
  snow on and off through the afternoon and evening, so hourly reports cannot separate the
  aircraft's contribution at the surface.

## 5. Calibration: what the test can and cannot tell apart

![Calibration](figures/calibration.png)

**Recall.** Synthetic streaks were added to the real radar of three quiet Denver days along
real flights' tracks, with physical timing: a supercooled layer 300-800 m deep at a height
the flight crossed, ice from each point reaching the radar beam after (its height above the
beam) / (a fall speed of 0.4-1.2 m/s), invisible for the first 10 minutes, 1.5-3 km wide
when it appears, and a wind error of 0.7 m/s. Of 426 such streaks, the detector found and
matched to some flight at p <= 0.05 (event recall) 16% at a peak of 5 dBZ, 52% at 10 dBZ,
60% at 15 dBZ and 63% at 20 dBZ; to the injected flight itself (attribution) 13%, 45%, 55%
and 57%. About one recovered streak in nine was credited to another aircraft, typically one
just ahead on the same path. Most misses are aircraft climbing or descending steeply
through the layer, whose ice reaches the beam piece by piece: 56% of those were recovered
against 71% of flights that crossed the layer nearly level, as the Denver seeders did at
their assigned approach altitudes. Census counts are therefore lower bounds, by a factor of
about 1.5-2 for streaks of 10-20 dBZ.

**Non-weather echo.** The first version of the test, using reflectivity alone, found
"aircraft streaks" on warm July days, when no aircraft can make snow (p = 0.005-0.01 on
three of the first six such radar-days), and one on a November day at Indianapolis. They
were echo that is not precipitation: ground clutter or turbines at a fixed spot along
Buffalo's southern arrival path, insects along summer convergence lines at Denver,
speckled clear-air echo at Indianapolis. The radars' correlation coefficient separates it
(snow and rain above ~0.97, such echo well below 0.9); masking echo below 0.90 removed
every one of those false alarms and left the Denver case intact.

**An independent review of the statistics** found that the null was too easy to beat:
neighbouring rotations are near-copies that match the same streaks, so a rotation scored
"as if real" was deflated by its own copies while the real flights, with no rotation within
20 deg, were not. Scoring every set only against sets at least 20 deg away, with a common
resolution floor, made the day's rank p-value calibrated in simulation (<= 0.05 on 4.0% of
null days). It also found that injection recovery counted matches to the wrong flight, and
several smaller bugs, all fixed (docs/methods.md lists them).

**How many independent rotations.** A day has only some 15-25 effectively independent
rotations. So no single day can be shown significant at the level a census-wide false
discovery rate needs, and on sparse days (a handful of streaks, rotations matching almost
nothing) even "beats all 200 rotations" is not rare: on 7 July 2025 at Detroit, four
arrivals matched one rain line along the final approach. Each day is therefore scored by
the margin z by which its flights beat their rotations, and its p-value is calibrated
against the same z of every rotation, pooled over all census days of similar chance-match
level. On 744 warm radar-days where no aircraft can make ice (every day of July and August
2025 and of July 2024 at the eight radars, flights below 13 000 ft in air above freezing)
this calibrated p-value was <= 0.01 on 7 days (7.4 expected) and <= 0.05 on 38 (37.2
expected).

**A corridor-preserving null cannot attribute streaks to aircraft.** Replacing the day's
flights with other days' flights at the same times of day, which fly the same corridors
through the same weather but cannot have seeded it, removes Denver along with the false
alarms (p = 0.10): on 7 of its 12 confirmed streaks (at the time), flights from 6-14 of the
28 other days fitted as well as the real seeder. Airline schedules repeat, so the same
flight numbers fly the same procedures at the same minutes every day; such a null keeps
clones of the seeders. Evidence that a streak is aircraft-made therefore rests on the
rotation test together with timing (echo-free air before the aircraft, echo after),
physics (supercooled air, snow-like fall speeds) and shape.

**Physical requirements** carry part of the weight: matches must be at or below -4 C, and
the ice must fall at snow speeds to reach the beam when it did; aircraft snow is light (at
most 35 dBZ). They removed one physically impossible Denver match (an aircraft at the
height of the beam itself).

**Time-reversed diagnostic.** Aircraft that pass 10-50 minutes after a scan cannot have
made its echo; at Denver they match no better than rotations (p = 0.12).

**Inspection.** Each candidate day's matched flights are also classified by fixed rules
(scripts/inspect_candidates.py): convincing when the streak is isolated (at most 30% of
the ring 2-8 km around it has echo), persists in two or more scans for the same flight,
and has a snow-like fall speed from below 16 000 ft; doubtful when it sits in widespread
echo (over 60% of the ring) or is physically odd. The rules were written with the first
reviewed days in view and checked against the Denver case; they are an aid, not a test.

## 6. Census

![Census results](figures/census.png)

**Scope.** Eight radars near large airports (Denver, Minneapolis, Chicago, Indianapolis,
Detroit, Cleveland, Buffalo, Pittsburgh), every day from November to March of the winters
2023-24, 2024-25 and 2025-26 plus March 2023: 3,736 radar-days. 848 were screened out
before testing, 474 with no echo, 366 with no cold, moist layer aloft and 8 with missing
radar data, which left 2,888 tested. A typical tested day had 915 flights below 30 000 ft
and 225 radar scans; 1.23 million streak objects were tested in all. Days whose 24-rotation
screen reached p <= 0.05 (355) were confirmed with 200 rotations. Nobody looked at any day
before its score was computed.

**How often.** If aircraft never made snow, 1% of tested days would reach a calibrated
p <= 0.01: 28.9 days. 58 did, an excess of 29 (90% interval 20-38). At p <= 0.05, 202 days
against 144.4 expected. Re-scored without streaks in families of parallel bands, the
census gives 64 and 204. The excess is not spread evenly: Chicago (15 days against 3.8
expected), Detroit (11 against 3.8) and Minneapolis (9 against 3.0) carry most of it;
Buffalo's 8 against 4.1 are mostly lake-effect band days that do not survive inspection,
and Cleveland, Indianapolis and Pittsburgh are at chance. By winter, 2023-24 shows no excess
(9 against 8.5), 2024-25 21 against 9.5 and 2025-26 27 against 9.5. In the main score no
single day reaches the bar a 5% false discovery rate over 2,888 days requires (section 5: a
day has only some 15-25 independent rotations); in the band-filtered score three do (below).
So the census establishes how often streaks occur, and the catalogue below which flights
made them, by different means.

**Supercooled cloud.** Aircraft can only make ice in cloud of supercooled droplets, and an
artefact of the method has no reason to care whether any is present. Split by the hours of
HRRR supercooled cloud water over the radar, days at p <= 0.05 are 3.4% with none (24 of
716), 7.9% with 1-3 hours (60 of 755), 6.1% with 4-8 hours (40 of 659) and 10.3% with 9-24
hours (78 of 758). After the first two winters showed this, a prediction about the colder
part of that covariate was written down before winter 2025-26 was analysed
(`docs/methods.md`): days without supercooled cloud at -10 C or colder would show no
excess, days with 4 or more hours of it would, and their share at p <= 0.05 would be the
larger. All three held: 23 of 398 days against 19.9 expected (no excess), 20 of 197 against
9.9 (one-sided Poisson p = 0.003), shares 5.8% against 10.2%.

**False alarms on warm days.** The 1% above is what a calibrated test should give; whether
this one does was measured. The same pipeline was run on every day of July and August 2025
and of July 2024 at the eight radars (744 radar-days), using only flight segments below
13 000 ft in air the HRRR puts above freezing, where no aircraft can make ice but traffic,
rain, insects and clutter are all present. The calibrated p-value was <= 0.01 on 7 of those
days (7.4 expected) and <= 0.05 on 38 (37.2 expected): the false-alarm rate is what it
claims to be. Measured against this rate instead of assumed, the winter excess at
p <= 0.01 is 29 days (90% interval 10-44) and at p <= 0.05 53 days (12-89); for the
band-filtered score, 35 (16-50) and 55 (14-91). The warm days' lowest p-value was 0.0009.
Three winter days in the band-filtered census have p = 0.00003, the smallest the
calibration can give, and pass a 5% false discovery rate on their own (q = 0.03): Chicago on
25 January and 2 February 2026 and Minneapolis on 14 February 2025. All three were judged
convincing by eye.

**What the streaks look like.** Every day at calibrated p <= 0.01 in either score (67) and
every day at p <= 0.05 with a flight the objective rules call convincing (34 more) was
reviewed by eye on a sheet of its eight best matches: 41 looked like aircraft streaks
(thin, in clear or light echo, along the drifted track, often growing over scans or tracing
the turn of an approach), 20 were ambiguous and 40 doubtful, mostly pieces of snow bands
lined up with the runways in use (`results/census/visual_review.csv`, with a note for each;
it also holds five days looked at during development, none convincing).

![Nine catalogue streaks](figures/census_gallery.png)

*Nine catalogue flights at six radars over three winters. Orange: the air each aircraft
flew through, carried by the HRRR wind to the scan time; dashed: its track over the previous
hour. Panel titles give the seeding altitude and temperature and the time since passage.*

**The catalogue.** On the 41 convincing radar-days, 316 flights matched a streak at
p <= 0.01, against 76 expected by chance from the rotations. Dropping those the objective
rules call doubtful (streak in widespread echo, among parallel bands, implausible fall speed
or seeding at 16 000 ft or above) leaves 163 flights on 35 radar-days (30 dates) at seven
radars: 16 in 2023-24, 99 in 2024-25 and 48 in 2025-26 (`data/census/seeders.csv`). They
are 155 jets, 4 turboprops, 3 piston aircraft and one of unknown type, much like the traffic
that crossed the seeding layer (96% jets). The commonest types are the 737-800 (19), CRJ-900
(14), E175 (13), CRJ-700 (12) and 737 MAX 8 (11); 24 are widebodies, from a 777-300ER to
MD-11 and 747-8 freighters. 114 were arriving, 48 departing and one in level flight.
They crossed the cloud at a median 7 300 ft (interquartile range 5 600-9 700 ft) and -13 C
(-17 to -11 C), and their ice implied fall speeds of 0.25-1.55 m/s (median 0.60). The
streaks were seen a median 33 minutes after passage, were 8-36 km long (median 13.5) and
about 2.3 km wide, and peaked at a median 1.0 mm/h of snow (water equivalent; at most 5).

**Arrivals, departures and temperature.** On catalogue days, 16.4% of jet arrivals that
entered echo-free air in the seeding layer left a streak along their drifted track (968 of
5,913), against 2.3% for the same flights rotated about the radar; for departures it was
2.5% (162 of 6,389) against 1.1%. Arrivals spend about twice as long in the layer, but a
Mantel-Haenszel odds ratio stratified by time in the layer, temperature and day still puts
them at 5.2 times the odds of departures (90% interval 4.1-6.5). What remains mixes aircraft
configuration (gear and flaps out, high lift at low speed) with airspeed, which this census
cannot separate. By temperature, arriving jets left streaks most often at -8 to -11 C
(25.6%, against 1.4% by chance), less often colder (12-18% below -11 C) and still well above
chance at -4 to -8 C (7.5% against 0.4%). Catalogue jets seeded at a median of -13 C, much
warmer than the -20 to -25 C often quoted for jets (`docs/background.md`). Days were selected
for having streaks, so these shares overstate the rate on an ordinary day; the comparisons
within them are the result. HRRR temperatures at the aircraft carry an error of 1-2 C.

**Without the visual review.** Taking every day at p <= 0.01 as a catalogue day instead
(no human judgement) gives 451 flights matched at p <= 0.01 against 89.5 by chance, and 213
catalogue flights on 47 radar-days at all eight radars, with the same picture: a median
seeding temperature of -13 C, 152 arrivals to 60 departures, fall speeds of 0.25-1.56 m/s.
The visual review changes which days are shown, not what they show.

## 7. Limitations

- **The chance model can be fooled.** Rotating flights about the radar breaks every link
  between where aircraft fly and where echo forms, including links that are not seeding
  (routing around storms, convergence lines near airports). On warm days such links give
  false alarms (7 of 744 at p <= 0.01, as many as a calibrated test should), so the
  census's false-alarm rate is measured on warm days at every radar, not assumed.
- **Corridors cannot be separated from aircraft by any placebo that keeps them.** Other
  days' flights fly the same procedures at the same minutes, so they fit real aircraft
  streaks almost as well as the aircraft that made them. Evidence that a streak is
  aircraft-made therefore rests on the rotation test, the timing (echo-free air before
  the aircraft, echo after), the physics (temperature, fall speed) and the shapes (turns of
  an approach pattern that weather does not draw).
- **Which aircraft.** In a stream on one path, the match is the best-fitting aircraft into
  echo-free air; an aircraft a few minutes ahead on the same path at the same height could
  have started the streak instead. Counts of seeding flights are conservative (followers
  through existing ice are not counted) and individual attributions within a stream are
  uncertain.
- **What can be seen.** Only seeding about 0.5-2 km above the lowest radar beam, close
  enough to fall into it within 10-50 minutes, and at least ~10 km of streak can be
  detected; aircraft descending steeply through a thin layer leave streaks too short to
  find. Rates are lower bounds.
- **Model fields.** HRRR often has no cloud water in thin supercooled stratus, and
  seeding temperatures are HRRR analysis temperatures at the aircraft's pressure level
  (typical error 1-2 C).
- **Snowfall.** Rates from reflectivity use one published Ze-S relation and are uncertain
  by a factor of a few, more for young, small crystals.
- **ADS-B coverage** from volunteer receivers changed over 2023-2026.
