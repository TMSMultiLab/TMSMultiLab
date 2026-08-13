# Constructing an open-source test dataset for MEP analyses

This directory will contain a series of scripts to retrieve and reformat open datasets for inclusion into a large, multi-lab open dataset.

The dataset will be used to test TMS, EMG, and MEP analysis code.

## Guidance

We don't want to store the data itself. We want to provide scripts that will allow users to download the data (semi-automatically), and recreate the full test dataset on their local machine.

One script per dataset.

For the first dataset:

* URL to download the data from
* Description of the dataset (plain text, context)
* Details of methods and analysis
* Data file (2D array, with one EMG epoch per row and one column per sample, e.g., 10000 rows and 2000 columns)
* Metadata file (2D array, with one row per row of the data file, and one column per metadata variable e.g., 10000 rows and 20 columns)
* Subjects file (2D array, with one row per subject and one column per subject metadata variable)

Future plans:
* Include all available metadata from the [TMS-RAT](https://tms-rat.org)  - see [notes here](https://github.com/TheHandLab/TMSMultiLab/edit/main/TestData/TMS-RAT.md)
* Include all available metadata from the NIBS-DAS
* Include all available metadata from NIBS-BIDS
* Convert all data to BIDS format

# TMS-RAT v1.0 metadata compliance

To ensure that analysed data are compliant with the TMS-RAT reporting recommendations, the following metadata should be included with all datasets (where relevant and available!).

Upper case letters A-N relate to the reporting categories in [Version 1.0](https://tms-rat.org/?page=rat&version=1)

A) Age, sex, handedness, height and weight

B) Medication, CNS, Neuro, Medical, Cognitive, Sleep, Posture

C) Studies, Population, Groups, Protocol, TimeofDay, SessionDuration/Interval, Attention, Control, Attenuation

D) Stimulator, CoilShape/Model, CoilDiameter, PulseWaveform

E) Neuronavigation, CoilLocationSkin, TargetCNS, CoilOrientation, CoilCurrentDirection, CoilStability

F) IntensityMethod, MachineIntensity, ParticipantIntensity

G) TrialInterval, NumberPulses, ResponseRaw

H) ThresholdLocation, ThresholdType, ThresholdCriterion, ThresholdAlgorithm, ThresholdValues

I) MotorSkills, PriorMotorActivity, MuscleActivity

J) EMGHardware, EMGSoftware, SkinPreparation, ElectrodeMaterial/Size, ElectrodeLocation, ElectrodeMontage, Sampling/Filtering, EMGMeasure, EMGAmplitudeNormalisation

K) AfferentStimulationHardware, AfferentElectrodes, AfferentLocation, AfferentTarget, AfferentDuration/Waveform, AfferentThresholdType, AfferentThreshold

L) ConditioningIntensityMachine, ConditioningIntensityParticipant, Condition-TestInterval(s), ConditionedUnconditionedOrder, ConditionedPulseNumber, ConditionedResponseRaw, ConditionedResponseRelative

M) SEPFrequency/Repetitions, EEGMontage, EEGFeature, EEGValues

N) BurstFrequency, WithinBurstPulses, BetweenBurstInterval, BurstNumber, BetweenTrainInterval, TrainNumber, TMSTaskInterval, RepeatNumber

# NIBS-DAS metadata compliance

# NIBS-BIDS metadata compliance
