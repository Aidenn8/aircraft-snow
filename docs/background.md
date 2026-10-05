# Background: how aircraft make it snow

This note summarises what was known before this project, with the sources checked. Numbers
are quoted from the papers' abstracts or the institutions' own summaries; where a figure
comes from a secondary source, that is said.

## The mechanism

Clouds of liquid droplets colder than 0 °C ("supercooled" clouds) are common in winter.
Pure droplets can stay liquid down to about −38 °C unless something starts freezing. An
aircraft can start it: air expanding over a propeller's blade tips or over a wing cools
by tens of degrees for a fraction of a second, enough to freeze droplets homogeneously.
The resulting ice crystals (aircraft-produced ice particles, APIPs) grow at the expense of
the surrounding droplets, because the vapour pressure over ice is lower than over water
(the Wegener-Bergeron-Findeisen process), and they fall out as snow. Seen from below, the
cloud gets a hole or a channel (hole-punch and canal clouds); seen by a weather radar, a
line of snow appears that follows the aircraft's path and drifts with the wind.

- Rangno & Hobbs (1983), *J. Appl. Meteor.* 22, 214: aircraft penetrations produced ice
  concentrations above 1000 per litre at −8 °C; the ice occupied a cylinder along the flight
  track about 300 m across after 5 minutes.
  [doi:10.1175/1520-0450(1983)022<0214:POIPIC>2.0.CO;2](https://doi.org/10.1175/1520-0450(1983)022<0214:POIPIC>2.0.CO;2)
- Woodley et al. (1991, 2003), *J. Appl. Meteor.*: controlled tests in supercooled fog and
  cloud between −5 and −12 °C. A Beech King Air "can produce APIPs readily" below −8 °C;
  some small piston aircraft did, others did not. Production is most likely at low
  temperature, high power, gear and flaps extended and low airspeed, i.e. in the approach
  and climb configurations. Mechanism: homogeneous freezing in air cooled around propeller
  tips. [2003 paper](https://doi.org/10.1175/1520-0450(2003)042<0640:AIPAAR>2.0.CO;2)
- Sassen (1991), *J. Appl. Meteor.* 30, 765: a King Air in a −28 to −31 °C altocumulus left
  242 ice crystals per litre along its track; the seeded volume spread about 1 m/s.
- Heymsfield et al. (2010), *BAMS* 91(6) ("Aircraft-induced hole punch and canal
  clouds: inadvertent cloud seeding"): two commercial turboprops documented producing holes
  and channels in altocumulus with NSF/NCAR C-130 radar, lidar and in-situ probes. Jets may
  do the same below their cruise levels. Commercial aircraft can affect clouds "at
  temperatures as much as 30 °C warmer than the −40 °C contrail formation threshold".
  [doi:10.1175/2009BAMS2905.1](https://doi.org/10.1175/2009BAMS2905.1)
- Heymsfield et al. (2011), *Science* ("Formation and spread of aircraft-induced
  holes in clouds"): model simulations show the ice growth drives vertical motions lasting an
  hour or more, which expand the holes; "regionally near major airports, additional
  precipitation can be induced". The accompanying NCAR release (30 June 2011) gives the
  thresholds usually quoted: turboprops at about −15 °C or colder, jets below about −20 to
  −25 °C; conditions that allow it occur up to about 6% of the time at the airports studied
  (more than 5% at Frankfurt, Paris CDG and Chicago O'Hare); one radar case showed an
  enhanced snowfall rate of close to an inch an hour; one Texas cloud deck had 92 holes, some
  lasting over four hours and up to 60 miles long.
  [doi:10.1126/science.1202851](https://doi.org/10.1126/science.1202851),
  [NCAR release](https://news.ucar.edu/4858/flights-can-inadvertently-generate-rain-and-snow-near-airports)
- Muraki, Rotunno & Morrison (2015), *J. Atmos. Sci.*: the expanding edge of a hole is a
  gravity-wave front.

## Observations of aircraft-made snow on radar

- Moisseev et al. (2019), *JGR Atmospheres* ("Inadvertent localized intensification of
  precipitation by aircraft"): eleven years of dual-polarisation radar near Helsinki-Vantaa.
  Localized areas a few to several tens of kilometres across where existing precipitation was
  6-14 times heavier than around them, linked to arriving and departing traffic; the snow in
  them had larger, more isotropic flakes. Proposed mechanism: aircraft ice boosting
  aggregation. The cases were found by inspection.
  [doi:10.1029/2018JD029449](https://doi.org/10.1029/2018JD029449)
- Mayes (2023), *Weather*: three snow stripes up to 10 km long near London Heathrow on
  23 January 2023, during freezing fog. [doi:10.1002/wea.4511](https://doi.org/10.1002/wea.4511)
- Iwasaki & Fujii (2025), *JGR Atmospheres* ("Formation of hole punch clouds"): 59 hole-punch
  clouds in Himawari-8 imagery matched to Flightradar24 tracks. Holes appear as lines along
  the tracks within 10-20 minutes of the aircraft passing, then widen at 1-2 m/s for over an
  hour; "not all aircraft that passed supercooled clouds generated" holes.
  [doi:10.1029/2025JD045417](https://doi.org/10.1029/2025JD045417)
- Denver, 25 January 2025 (local evening): a radar pattern of snow near Denver International
  Airport between about 6:10 and 7 pm MST, linked in news coverage to United flight 5528
  (operated by SkyWest) and about a dozen other arrivals using a north-south runway.
  [Interesting Engineering, 28 Jan 2025](https://interestingengineering.com/culture/skywest-flight-triggers-rare-snowfall),
  citing the Washington Post and View from the Wing.

## Related: snow from industry

Toll et al. (2024), *Science*: satellite evidence of supercooled clouds glaciating
downwind of industrial aerosol hot spots at −10 to −24 °C, with daily snowfall accumulations
up to 15 mm. Industrial snow is a possible false positive for this project, but it starts at
fixed stacks rather than along flight paths. [doi:10.1126/science.adl0303](https://doi.org/10.1126/science.adl0303)

## What has not been done

A search of citing literature (OpenAlex, October 2026) for all of the above found case
studies and one single-airport inspection (Helsinki), but no automated detector and no
census: nobody has measured, across many airports and winters, how often aircraft make
radar-visible snow, which aircraft do it, at what temperatures, and how much snow results.
The temperature thresholds in use (turboprops about −10 to −15 °C, jets about −20 to
−25 °C) come from a handful of cases, field tests with research aircraft, and sounding
climatologies, not from counting real events.

## What the open data allow

- NEXRAD WSR-88D radars: Level II and Level III files on AWS (Unidata/NOAA buckets),
  minutes apart, 250 m range gates at 0.5 degree azimuth spacing at the lowest tilt.
- adsb.lol: open (ODbL) global ADS-B history since 10 March 2023, a position per aircraft
  every 10-20 s with pressure altitude; aircraft type from the tar1090 database.
- HRRR analyses: hourly 3 km temperature, humidity, cloud water and wind profiles
  (University of Utah Zarr archive on AWS).
- METARs and pilot reports (Iowa Environmental Mesonet).
- Radiosondes: NOAA IGRA2. The Denver station's archive ends on 8 July 2022 in both IGRA2
  and the Iowa Environmental Mesonet, so the Denver case relies on HRRR for profiles.
