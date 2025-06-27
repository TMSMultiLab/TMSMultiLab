# Code to analyse TMS data

## EMG_filter.m
<code>MEP.m</code> accepts an m x n array of time-series data, where individual channels of data are in columns. needs at minimum the sampling frequency (to perform a 2nd order, zero phase 20-500Hz bandpass). it returns the same data after filtering. options include order (default=2), lower frequency cutoff (default=20), upper frequency cutoff (default=500).

<code>filtered_data = EMG_filter(data,samplehz,order,lowhz,highhz);</code>

within a script, specify, eg:

<code>samplehz=4000;</code> - sampling frequency in samples per second (at least 2000 is recommended)

<code>order=1;</code> - order of filter, which is used twice, so 1 -> 2nd order, 2 -> 4th order (default = 1, 2nd order)

<code>lowhz=20;</code> - lower cut-off of the filter (default 20)

<code>highhz=500;</code> - upper cut-off of the filter (default 500)


## MEP.m
<code>MEP.m</code> accepts a single column of time-series data and needs the sampling frequency and the time that the stimulation was given. options can also be given. it returns a structure of measurements.

example:

<code>data=randn(8000,1);</code> - single column of data (1000ms before and after the pulse is recommended)

<code>samplehz=4000;</code> - sampling frequency in samples per second (at least 2000 is recommended)

<code>pulsetime=1000;</code> - time in ms that the stimulus was presented (1000ms is recommended)

<code>[mep,options]=MEP(data,samplehz,pulsetime);</code>

the <code>mep</code> structure contains fields including the peak-to-peak amplitude, <code>mep.amp</code>
