# BioMetrics data relevant to TMS studies
The data and code in this folder provide estimates of important variables used in TMS research, for example the latencies and amplitudes of evoked potentials (sensory or motor), and the relationship between heights, ages, sex, and head sizes of participants. Optimal use of this information is likely to decrease between-participant variation in basic TMS methods.

## HandLab Data
Data preceded by 'HandLab' were collected starting in 2011 in the HandLab. For each participant in a TMS study, we measure the head from nasion to inion and between the pre-auricular points.

More recently we also measure the head circumference (nasion - pre-auricular point - inion), and we are also starting to systematically measure participant height, weight, and arm length.

These biometric data are likely useful in
* training researchers to measure heads
* finding typical locations on the scalp
* correlating MEP parameters such as amplitude and latency with height and arm length

Files:
1) <b>HandLab_TMSHeads.csv</b>  - the data, in comma-separated values:
   
   headid - unique ID number for this row
   
   headtype - mostly 'head', though some 'arm' measurements too
   
   participantid - unique person identifier. rows with the same participantid come from the same human
   
   headdate - date and time that the head was measured
   
   nasioninion - distance in cm between nasion and inion, measured through the vertex
   
   intertragal - distance in cm between the two pre-auricular points, measured through the vertex
   
   nasionearinion - distance in cm beween the nasion and inion, passing just above the pre-auricular point (above the ear, but as low as possible, behind the superior pinna)
   
   armlength - length in cm between the tip of the most distal outstretched finger and the acromion (https://en.wikipedia.org/wiki/Acromioclavicular_joint)
   
   wristcirc - during peripheral nerve experiments we sometimes record the circumference of the wrist, in cm

2) <b>HandLab_TMSHeads.m</b> - Matlab script (will also work with minor changes using Octave (https://octave.org/) to load and parse the data and produce the graphs shown on the Head measurement page (https://github.com/TMSMultiLab/TMSMultiLab/wiki/Head-measurement)

3) <b>HandLab_TMSSites.m</b> - Matlab script to plot the locations of muscles targeted onto a map of the scalp, centred on the vertex, Cz, as shown on the [primary motor cortex](https://github.com/TMSMultiLab/TMSMultiLab/wiki/Primary-motor-cortex) page.
