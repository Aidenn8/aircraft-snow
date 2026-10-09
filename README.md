# aircraft-snow

When an airliner climbs or descends through a cloud of supercooled droplets, air cooling
over its wings and propellers can freeze some of them. The ice grows, falls out as snow and
drifts with the wind. On weather radar the result is a thin line of snow echo that traces
the aircraft's path, offset downwind.

The effect has been documented in case studies (hole-punch and canal clouds, single radar
events near Denver, Helsinki and London). This repository is the first attempt at a census:
an automated detector that finds these streaks in NEXRAD radar data, matches each one to the
aircraft that made it using open ADS-B flight tracks, and measures how often it happens, at
what temperatures, and which aircraft do it, with false-alarm rates measured rather than
assumed.

**Results** (`REPORT.md`): over 2,888 tested winter radar-days at eight radars (November
2023 to March 2026, plus March 2023), 58 days show flight-streak matches beyond chance at
p <= 0.01 where 29 would be false alarms, a rate measured on 744 warm days on which no
aircraft can make ice. The excess follows supercooled cloud over the radar, as predicted
before the last winter was analysed. A catalogue of 163 flights on days judged convincing
by eye (seven radars; `results/census/seeders.csv`) shows seeding
at a median of 7 300 ft and -13 C, and arriving jets leaving streaks at about five times
the odds of departing ones. The Denver case of 25 January 2025 is worked through in detail:
16 flights matched streaks along their wind-carried tracks against 3.6 expected by chance.

## How it works

1. **Flights** (`snow/adsb.py`, `snow/flights.py`): every aircraft position over North
   America below 30 000 ft from the open adsb.lol archive (one position per aircraft per
   10 s), split into flights, identified with the tar1090 aircraft database.
2. **Radar** (`snow/radar.py`, `snow/event.py`): NEXRAD Level III lowest-tilt reflectivity
   from AWS, resampled to a 500 m grid around the radar. Aircraft echoes are removed where
   the beam saw each aircraft.
3. **Wind** (`snow/hrrr.py`, `snow/epoch.py`): hourly HRRR analysis profiles; the drift of
   seeded air is the mean wind between the radar beam and the aircraft's level.
4. **Streaks** (`snow/streaks.py`): thin, long echo objects found by morphology alone, with
   fixed clutter removed.
5. **Matching** (`snow/match.py`): each streak is compared with every flight that passed
   10-50 minutes earlier, its track carried by the wind. A match needs the streak to lie along
   the drifted track (coverage and precision), the air to have been echo-free before the
   aircraft arrived, and physically possible seeding: air at or below -4 C, and a height above
   the radar beam that ice falling at snow speeds (0.2-2 m/s) would cover in the observed time.
   Chance is measured by repeating everything with all flights rotated about the radar.
6. **Day test** (`snow/match.py: day_statistic`): the evidence from all matched flights of a
   radar-day, against the same statistic for every rotated flight set; across the census,
   days are selected at a 5% false discovery rate after calibrating the statistic's tail on
   the rotations themselves (`scripts/census_summary.py`).
7. **Calibration**: synthetic streaks injected into real radar data (`snow/inject.py`,
   `scripts/injection_recovery.py`) measure recall; warm-season days where seeding is
   impossible, translated instead of rotated flights (`scripts/shift_null.py`) and a
   time-reversed test measure false alarms.

A second test (`snow/event.py`), which follows every track point's air parcel and compares
the change in reflectivity on the track with lines beside it, turned out not to be
calibrated on winter days and is not used in the census (see `docs/methods.md`).

See `docs/background.md` for the science and prior work and `docs/methods.md` for the method
and the lessons learned building it.

## Running it

```bash
python -m venv .venv && .venv/bin/pip install -e '.[test]'
.venv/bin/python scripts/fetch_reference.py          # aircraft database
.venv/bin/python scripts/build_sites.py              # airports and radars
.venv/bin/python -m snow.adsb 2025-01-25:2025-01-26  # ADS-B days (~25-60 s each)
.venv/bin/python scripts/run_days.py KFTG 2025-01-25 2025-01-25
.venv/bin/python scripts/review_day.py KFTG 2025-01-25
.venv/bin/python -m pytest
```

Long runs on a laptop without a fan go through `tools/thermal/guarded_run.py`, which pauses
the job when the machine gets hot (start `tools/thermal/daemon.sh` first).

## Data and licences

- ADS-B: [adsb.lol](https://www.adsb.lol/) historical data, Open Database License (ODbL).
- Radar: NOAA NEXRAD Level II/III via the AWS Open Data program (Unidata buckets).
- Weather model: NOAA HRRR, University of Utah HRRR-Zarr archive on AWS.
- Aircraft types: [tar1090-db](https://github.com/wiedehopf/tar1090-db).
- Surface observations: Iowa Environmental Mesonet.
