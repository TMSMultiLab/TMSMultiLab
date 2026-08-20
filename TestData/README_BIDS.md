%% HOW TO BIDSIFY YOUR DATASET__________________________________
%
% GUIDANCE
% https://bids-specification.readthedocs.io/en/stable/modality-specific-files/electromyography.html
% Do not include white spaces in file names: They make scripting harder.
% Use only letters, numbers, hyphens, and underscores: Some operating systems cannot handle special characters.
% Do not rely on letter case (UPPERCASE and lowercase): For some operating systems a is the same as A.
% Use separators and case in a systematic and meaningful way: thisIsCamelCase, this_is_snake_case
% key1 - value1 _ key2 - value2 _ suffix .extension
% Suffixes are preceded by an underscore
% Entities are composed of key-value pairs separated by underscores
% There is a limited set of suffixes for each data type (anat, func, eeg, ...)
% For a given suffix, some entities are required and some others are [optional].
% Entity keys and suffixes can only contain letters and/or numbers.
% Entity values can only contain letters and/or numbers and/or the "+" character.
% Entity key-value pairs have a specific order in which they must appear in filename.
% Some entities key-value can only be used for derivative data
%
%
% VALIDATION
% check .json file formatting by opening in a web browser - should appear as a nice table without any errors


%% ALL BIDS REQUIRE THE FOLLOWING FILES__________________________
%

%% DATASET_DESCRIPTION FILES_____________________________________
% https://bids-specification.readthedocs.io/en/stable/modality-agnostic-files/dataset-description.html
%
%	dataset_description.json
%		Name				REQUIRED
%		BIDSVersion			REQUIRED
%
%		HEDVersion			RECOMMENDED
%		DatasetType			RECOMMENDED
%		License				RECOMMENDED
%		Authors				RECOMMENDED
%		GeneratedBy			RECOMMENDED	[object]		provenance of the dataset, eg derivitives generated manually or by code
%		SourceDatasets			RECOMMENDED	[objects]		URLs  / DOIs to the original data
%
%		DatasetLinks			OPTIONAL				BIDS URIs
%		Keywords			OPTIONAL				free text
%		Acknowledgements		OPTIONAL				other people/institutions not included in Authors or funding
%		HowToAcknowledge		OPTIONAL				how should re-users of the dataset acknowledge the authors
%		Funding				OPTIONAL				list of funding sources - grant numbers
%		EthicsApprovals			OPTIONAL				ethics approval committee / ID; protocol numbers
%		ReferencesAndLinks		OPTIONAL				references to additional information about the data / URLs
%		DatasetDOI			OPTIONAL				full URI to the doi (https://doi.org/...)

% 	README[.md|.rst|.txt]			REQUIRED				plain text description of the dataset
%	CITATION.cff				OPTIONAL
%	CHANGES					OPTIONAL				must follow CPAN changelog conventions
%	LICENSE[.md|.rst|.txt]			OPTIONAL				plain text descriptipon; must match the dataset_description licence field


%% DATA SUMMARY FILES____________________________________________
% https://bids-specification.readthedocs.io/en/stable/modality-agnostic-files/data-summary-files.html
%	participants.tsv			RECOMMENDED				eg age, sex, handedness, species and strain
%		participant_id			REQUIRED				first column, sub-01..., unique single row per P
%
%		species				RECOMMENDED	homo sapiens		
%		age				RECOMMENDED				number in years. for ages>88, use 89+; avoid months/weeks unless eg babies/animals
%		sex				RECOMMENDED				"male" "female" "other"
%		handedness			RECOMMENDED				"L", "R", "A" (case-insensitive, also 'left' 'right' 'ambidextrous'
%		strain				RECOMMENDED
%		strain_rrid			RECOMMENDED
%
%		HED				OPTIONAL
%		[additional columns]		OPTIONAL				any other columns must be described in meta-data file
%
%	participants.json			RECOMMENDED				for each column above, free text description, eg:
%		{
%			"age": {
%			"Description": "age of the participant",
%			"Units": "year"
%		}
%  
% Samples file					OPTIONAL				for biological samples
%	samples.tsv
%	samples.json
%
% Scans file					OPTIONAL				one row per datafile
%	*_scans.tsv
%		filename			REQUIRED				relative paths, first column
%		acq_time			OPTIONAL				YYYY-MM-DDThh:mm:ss[.000000][Z|+hh:mm|-hh:mm]; YEAR can be shifted to <=1925 to preserve anonymity
%		HED				OPTIONAL				
%
%       *_scans.json
%
% Sessions file
%	*_sessions.tsv				OPTIONAL				if multiple sessions per participant
%		session_id			REQUIRED				1st column, unique ID/number for this session, one file per participant
%
%		pathology			RECOMMENDED
%
%		acq_time			OPTIONAL
%		HED		EMGCoordinateSystem		OPTIONAL
%		[additional columns]		OPTIONAL				any other columns must be described in meta-data file


% PHENOTYPE AND ASSESSMENT FILES___________________________				eg survey results / quesitonnaires too big for participants files
% https://bids-specification.readthedocs.io/en/stable/modality-agnostic-files/phenotypic-and-assessment-data.html
%	phenotype/<measurement_tool_name>.tsv						
%	phenotype/<measurement_tool_name>.json


% CODE_____________________________________________________
% https://bids-specification.readthedocs.io/en/stable/modality-agnostic-files/code.html
%
%	code/*					OPTIONAL		


% EVENTS___________________________________________________
% https://bids-specification.readthedocs.io/en/stable/modality-agnostic-files/events.html
% 
%       [modality]/<matches>_events.json	OPTIONAL				see individual modalities for guidance
%       [modality]/<matches>_events.tsv


%% ------------------------------------------------------------------------------------------------- %%

%% EMG BIDS (BEP 042) REQUIRES THE FOLLOWING FILES_______________

%% EMG FILE NAMING_________________________________________________________
% sub-<label>[_ses-<label>]_task-<label>[_acq-<label>][_run-<index>][_recording-<label>]_
% eg:
% sub-01_ses-01_task-rest_acq-M125_run2_recording-highres.edf
% sub-01_ses-01_task-pegboard_acq-SMG45.edf
%
% _acq - main parameters?
% _run-1, _run-2 etc - repeats of the same acquisition parameters and conditions
% _recording-sEMG, _recording-HDEMG - when multple recordings taken at the same time


% https://bids-specification.readthedocs.io/en/stable/modality-specific-files/electromyography.html
%
%	*_emg.edf
%
%
%	*_emg.json
%		EMGPlacementScheme		REQUIRED 	"Other"			Must be one of: "ChannelSpecific", "Measured", "Other".
%		EMGPlacementSchemeDescription	REQUIRED 	"[description]"		if 'Other' above - description of general scheme, not muscle-specific
%		EMGReference			REQUIRED 	"ChannelSpecific"	(Bipolar is for fixed width electrodes
%		SamplingFrequency		REQUIRED	4000
%		PowerLineFrequency		REQUIRED	50
%		RecordingType			REQUIRED	"epoched"		Must be one of: "continuous", "epoched", "discontinuous".
%		SoftwareFilters			REQUIRED	"n/a"			"n/a" or object with details
%		TaskName			REQUIRED	"rest" / "pegboard"	eg "faces n-back", "head nodding", "treatment", "control", "sleep"
%		[additional columns]		OPTIONAL				any other columns must be described in meta-data file
%
%		TaskDescription			RECOMMENDED
%		Instructions			RECOMMENDED				eg "seated", "reclined"
%		EMGChannelCount			RECOMMENDED	4			integer
%		HardwareFilters			RECOMMENDED	{"Highpass filter": {"Cutoff (Hz)": 10"}, "Lowpass filter": {"Cutoff (Hz)": 2000"}}.
%		RecordingDuration		RECOMMENDED				integer
%		EpochLength			RECOMMENDED	2			duration in seconds
%		Manufacturer			RECOMMENDED	"AD Instruments"
%		ManufacturersModelName		RECOMMENDED	"PowerLab 16/30"
%		SoftwareVersions		RECOMMENDED	"LabChart 8"
%		DeviceSerialNumber		RECOMMENDED
%		ElectrodeManufacturer		RECOMMENDED	
%		ElectrodeManufacturersModelName	RECOMMENDED	
%		InstitutionName			RECOMMENDED	"University of Birmingham"
%		InstitutionAddress		RECOMMENDED	"Edgbaston, Birmingham, B15 2TT, UK"
%		InstitutionalDepartmentName	RECOMMENDED	"School of Sport, Exercise and Rehabilitation Sciences"
%
%		ElectrodeMaterial		OPTIONAL	"Ag/AgCl"		eg Tin, Ag/AgCl, Gold
%		ElectrodeType			OPTIONAL	"nipple"		eg cup, ring, clip-on, wire, needle
%		EMGGround			OPTIONAL				eg "right radial styloid process"
%		Gain				OPTIONAL				gain between electrode and the DAQ
%		InterelectrodeDistance		OPTIONAL	2.5			distance (in cm?) between electrode centres
%		Preamplification		OPTIONAL				pre-amp within the sensor?
%		SkinPreparation			OPTIONAL	"alcohol wipe" 		eg "alcohol wipe" or "abrasive gel"
%		SubjectArtefactDescription	OPTIONAL	"n/a"			eg "Vagus Nerve Stimulator", "non-removable implant"; "n/a" - no major artefacts
%		TriggerChannelCount		OPTIONAL	2			N channels with digital or analogue triggers
%
%
%	*_channels.tsv
%		name				REQUIRED	"time" / "right FDI"	1st column, unique name for this channel
%		type				REQUIRED	LATENCY	/ EMG / TRIG	2nd column, UPPER CASE, MUST BE: ECG / EMG / EOG / HEOG / LATENCY / MISC / REF / SYSCLOCK / TRIG / VEOG
%		units				REQUIRED	"s" / "V"		3rd column, physical units
%
%		signal_electrode		RECOMMENDED				for bipolar devices - same as channel name?
%		reference			RECOMMENDED	"bipolar"		
%		group				RECOMMENDED	"device1"		eg when electrodes are clustered together within a device
%		target_muscle			RECOMMENDED	"left fdi"; "flexors"	
%		placement_scheme		RECOMMENDED	"other"			must be: "measured", "other"
%		placement_description		RECOMMENDED				description
%		interelectrode_distance		RECOMMENDED	
%
%		description			OPTIONAL				free text, eg: n/a, stimulus, response, skin conductance, battery status
%		sampling_frequency		OPTIONAL	4000			(not needed if specified in _emg.json)
%		low_cutoff			OPTIONAL	10			number
%		high_cutoff			OPTIONAL	2000			number
%		notch				OPTIONAL	"n/a"			eg 50, 60, [60, 120, 180]
%		status				OPTIONAL	"good"			must be "good" or "bad"
%		status_description		OPTIONAL				free text for why it's bad
%		[additional columns]		OPTIONAL				any other columns must be described in meta-data file
%
%
%	*_electrodes.tsv
%		name				REQUIRED
%		x				REQUIRED
%		y				REQUIRED
%
%		coordinate_system		RECOMMENDED
%		type				RECOMMENDED
%		material			RECOMMENDED
%		impedance			RECOMMENDED
%		group				RECOMMENDED	"device1"		name and group combination must be unique; matches *_channels.tsv
%
%		z				OPTIONAL
%		[additional columns]		OPTIONAL				any other columns must be described in meta-data file
%
%
%	*_coordsystem.json								REQUIRED if *_electrodes.json used
%		EMGCoordinateSystem		REQUIRED	"other"			"Other" is the only valid value until standardized names become available for EMG coordinate systems
%		EMGCoordinateUnits		REQUIRED	"n/a"			Must be one of: "m", "mm", "cm", "percent", "n/a".
%		EMGCoordinateSystemDescription	REQUIRED	[description]		
%
%		ParentCoordinateSystem		OPTIONAL
%		AnchorCoordinates		OPTIONAL				(required if parent system used)
%		AnchorElectrode			OPTIONAL				(required if parent system used)
%
%
%	*_photo.<extension>			OPTIONAL				photo of the montage