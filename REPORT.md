# Aircraft-made snow: an automated census from open radar and flight data

*Draft. The Denver case study and the method are complete; the multi-airport census and its
calibration are running, and their sections will be filled in from the results.*

## Summary

Aircraft climbing or descending through clouds of supercooled droplets can freeze some of
them; the ice grows, falls and drifts with the wind, and weather radar sees a narrow band
of snow that traces the flight path, displaced downwind. This had been shown in individual
cases. This project builds an automated detector for these streaks in NEXRAD radar data,
matches each streak to the single aircraft that made it using open ADS-B tracks, and
measures how often chance alone would produce such a match.

On the evening of 25 January 2025 at Denver International Airport, the detector finds
streaks from 14 identified jets at a false discovery rate of about 5%, the day's combined
evidence exceeding all 200 placebo flight sets (gamma-tail p = 2e-4). Of the jets whose
first passage through the -17 to -20 C cloud layer could be followed, 62% left a streak,
against 10% expected by chance.

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
| NEXRAD Level III N0B (AWS) | lowest-tilt reflectivity, 0.5 deg x 250 m |
| NEXRAD Level II (AWS) | dual-polarisation moments for the case study |
| HRRR analyses (Utah Zarr archive) | hourly temperature, humidity, wind profiles |
| tar1090-db | aircraft type and engine class |
| Iowa Environmental Mesonet | METARs |

## 3. Method

See `docs/methods.md` for the full description and the lessons from building it. In brief:

1. Streak-shaped echo objects (long, thin, at most 35 dBZ) are found in each radar scan by
   morphology alone.
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

![Denver radar snapshots](figures/denver_case.png)

- **Conditions.** Overcast at 4 800 ft above the airport (~3.1 km above sea level), -8 C at
  the surface, light easterly upslope wind; HRRR puts the cloud layer at about -16 to -19 C
  under a warm inversion at 3.6 km.
- **Streaks.** From about 23:15 to 03:30 UTC the lowest radar tilt shows sharp bands of
  snow echo, 9-30 km long and 2-3 km wide, tracing the downwind, base and final legs of the
  arrival pattern for the north-flow runways, and other paths, displaced downwind.
- **Aircraft.** At p <= 0.005, 14 flights: CRJ-200, E175, A319, A320, 737-800, 737 MAX 8,
  737-900, 757, 787-9 and a Citation Latitude; mostly arrivals, seeding at 6 800-13 200 ft
  and -12 to -18.5 C (HRRR), streaks visible 13-46 minutes after passage.
- **Drift.** The drift that best aligns tracks with echo agrees with the HRRR wind to about
  1 m/s.
- **Particles.** At the same reflectivity, streak snow has lower differential reflectivity
  than the surrounding light snow (median 0.22 vs 0.44 dB), as Moisseev et al. (2019) found:
  aggregated, more isotropic flakes.
- **How many jets did it.** Of 48 jets whose first passage through the -17 to -20 C layer
  could be followed, 30 (62%) left a streak along their drifted track; the same flights
  rotated about the radar give 10%.
- **Ground.** Nearby stations reported light snow on and off all evening, so hourly METARs
  cannot isolate the aircraft contribution at the surface.

## 5. Calibration

*To be completed from the July null runs, the January runs, the time-reversed diagnostic
and the injection-recovery tests.*

## 6. Census

*To be completed.*

## 7. Limitations

- Detection is limited to seeding within roughly 1-1.5 km above the lowest radar beam:
  higher seeding takes longer than the 50-minute matching window to fall into view.
- HRRR often has no cloud water in thin supercooled stratus; seeding temperatures are HRRR
  analysis temperatures at the aircraft's pressure level (typical error 1-2 C).
- In a stream of aircraft on one path, only the first through new cloud can be matched; the
  census counts seeding flights conservatively.
- Snowfall rates from reflectivity use one published Ze-S relation and are uncertain by a
  factor of a few.
- ADS-B coverage from volunteer receivers changed over 2023-2026.
