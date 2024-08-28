The data in this folder were collected starting in 2011 in the HandLab. For each participant in a TMS study, we measure the head from nasion to inion and between the pre-auricular points.
More recently we also measure the head circumference (nasion - pre-auricular point - inion), and we are also starting to systematically measure participant height, weight, and arm length.
These biometric data are likely useful in a) training researchers to measure heads; b) finding typical locations on the scalp; c) correlating MEP parameters such as amplitude and latency with height and arm length.

Files:
1) HandLab_TMSHeads.csv  - the data, in comma-separated values:
   headid - unique ID number for this row
   headtype - mostly 'head', though some 'arm' measurements too
   participantid - unique person identifier. rows with the same participantid come from the same human
   headdate - date and time that the head was measured
   nasioninion - distance in cm between nasion and inion, measured through the vertex
   intertragal - distance in cm between the two pre-auricular points, measured through the vertex
   nasionearinion - distance in cm beween the nasion and inion, passing just above the pre-auricular point (above the ear, but as low as possible, behind the superior pinna)
   armlength - length in cm between the tip of the most distal outstretched finger and the acromion (https://en.wikipedia.org/wiki/Acromioclavicular_joint)
   wristcirc - during peripheral nerve experiments we sometimes record the circumference of the wrist, in cm

2) HandLab_TMSHeads.m - Matlab script (will also work with minor changes using Octave (https://octave.org/) to load and parse the data and produce the graphs shown on the Head measurement page (https://github.com/TMSMultiLab/TMSMultiLab/wiki/Head-measurement)
3) 
