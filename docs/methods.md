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
flights that do not cross cloud. This test is kept as an independent, assumption-light
check, but is not the main detector.

## Test 2: streak objects matched to single flights (streaks.py, match.py)

1. Streak objects are found in each scan without reference to aircraft: Gaussian smoothing
   (0.75 km), white top-hat with a 4 km disk (removes broad precipitation), threshold
   (top-hat >= 6 dB and >= 3 dBZ), cleanup, then shape by skeleton: main path >= 8 km,
   mean width <= 3.5 km, branchiness <= 2.5. Skeleton measures accept the U-shaped
   streaks left by approach patterns.
2. Fixed echoes are removed: objects lying mostly on pixels that are inside streak objects
   in >= 40% of the day's scans (ground clutter on terrain edges).
3. Each object is compared with every flight that passed 10-50 min earlier. Tracks are
   interpolated to 250 m spacing and drifted with the HRRR drift plus a perturbation of
   up to +-1 m/s (0.5 m/s steps). Coverage (share of the object's centre line within 1 km
   of the drifted track) and precision (share of the drifted track near the object that is
   on it) are combined as F1.
4. Newness: the object, moved back with the drift, must have been echo-free (>= 70% of
   samples < 0 dBZ) in the scans 3-30 min before the flight passed. In a stream of aircraft
   the match is therefore the first aircraft through new cloud.
5. Chance: the whole procedure is repeated with all flights rotated about the radar
   (24 rotations to screen a day, 200 to confirm). A day is scored by E = sum over objects of
   -log10 p; its null distribution treats each rotation in turn as the real set.

Denver (full day): 896 streak objects, 81 matched; E = 104 against a null median of 43 and
95th percentile of 70; no rotation reached the real value.

## Calibration

- Synthetic streaks injected along real flights' drifted tracks into real radar on quiet
  days (inject.py: Gaussian cross-section 0.5 km growing 0.2 m/s, onset 15 min, peak 35 min,
  gone at 90 min, wind error 0.7 m/s per component). Recall on three quiet Denver days
  (120 trials): 23% at 5 dBZ peak, 93% at 10 dBZ, 97% at 15 dBZ, 79% at 20 dBZ.
- False alarms: the same pipeline on July days restricted to tracks below 13 000 ft, where
  the air is above freezing and aircraft cannot seed (null mode).

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
