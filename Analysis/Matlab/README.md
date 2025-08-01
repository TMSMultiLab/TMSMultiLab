# Code to analyse TMS data

## spike_artefact.m
<code>spike_artefact.m</code> attempts to remove the electrical or magnetic artefacts associated with nerve or brain stimulation. it accepts an mx1 array of m samples of data and requires the sampling frequency and the stimulation time. uses a double-exponential fitting function to esimate the shape of the spike artefact over the period between 2ms and 75ms after the stimulus, and subtracts this function from the data. before subtraction, the exponential function is multiplied by a square-root function to progressively decrease the influence of the exponentials (and avoid step artefacts at the end). after subtraction, the data are linearly interpolated between 0 and 2ms (the 'non-recoverable period' of Erez et al., 2010). various options are available (need to put these all in <code>opts</code>). the function outputs a structure containing the cleaned (and filtered) data and (optionally) plots data to illustrate the artefact-removal process.

<code>[spike, fitresult, gof] = spike_artefact(data, samplehz, stimtime [,artefactwindow] [,recovery] [,scale] [,exclude] [,plot_artefact]);</code>

within a script, specify, eg:

<code>data = randn(8000,1);</code> - single column of data

<code>samplehz = 4000;</code> - sampling frequency in samples per second (at least 2000 is recommended)

<code>stimtime = 0.2;</code> - time in seconds that the stimulus was presented (at least 100ms of pre-stimulus baseline is recommended)

requires:
* EMG_filter (TMSMultiLab)
* fit (Matlab)

used by:
* cutaneomotor (TMSMultiLab)


## EMG_filter.m
<code>EMG_filter.m</code> accepts an mxn array of time-series data, where individual channels of data are in n columns. needs at minimum the sampling frequency (to perform a 2nd order, zero phase 20-500Hz bandpass). it returns the same data after filtering. options include order (default=2), lower frequency cutoff (default=20), upper frequency cutoff (default=500).

<code>filtered_data = EMG_filter(data,samplehz,order,lowhz,highhz);</code>

within a script, specify, eg:

<code>data = randn(8000,4);</code> - four columns of data

<code>samplehz = 4000;</code> - sampling frequency in samples per second (at least 2000 is recommended)

<code>order = 1;</code> - order of filter, which is used twice, so 1 -> 2nd order, 2 -> 4th order (default = 1, 2nd order)

<code>lowhz = 20;</code> - lower cut-off of the filter (default 20)

<code>highhz = 500;</code> - upper cut-off of the filter (default 500)


## cutaneomotor.m
<code>cutaneomotor.m</code> takes an mxn array of raw EMG data (m samples x n repeitions) and estimates changes in rectified mean EMG signal following a stimulus. designed to analyse the cutaneo-muscular reflex, it could be generalised to silent periods, h-reflexes or MEPs.


## MEP.m
<code>MEP.m</code> accepts a single column of time-series data and needs the sampling frequency and the time that the stimulation was given. options can also be given. it returns a structure of measurements.

example:

<code>data = randn(8000,1);</code> - single column of data (1000ms before and after the pulse is recommended)

<code>samplehz = 4000;</code> - sampling frequency in samples per second (at least 2000 is recommended)

<code>stimtime = 1000;</code> - time in ms that the stimulus was presented (1000ms is recommended)

<code>[mep, options] = MEP(data, samplehz, pulsetime);</code>

the <code>mep</code> structure contains fields including the peak-to-peak amplitude, <code>mep.amp</code>
