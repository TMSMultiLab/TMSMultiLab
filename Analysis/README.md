# code to analyse TMS data

## MEP.m
<code>MEP.m</code> accepts a single column of time-series data and needs the sampling frequency and the time that the stimulation was given. options can also be given. it returns a structure of measurements.

example:

<code>data=randn(8000,1);</code> - single column of data (1000ms before and after the pulse is recommended)

<code>samplehz=4000;</code> - sampling frequency in samples per second (at least 2000 is recommended)

<code>pulsetime=1000;</code> - time in ms that the stimulus was presented (1000ms is recommended)

<code>[mep,options]=MEP(data,samplehz,pulsetime);</code>

the <code>mep</code> structure contains fields including the peak-to-peak amplitude, <code>mep.amp</code>
