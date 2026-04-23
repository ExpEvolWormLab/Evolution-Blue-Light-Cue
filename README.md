# CueReliability
repository to archive raw data, modeling scripts and results of the paper on cue reliability Teotónio and Proulx 2025

FitnessAnalysis.Rmd contains modeling scripts and results, including script for the figures

cue_data contains the assay fitness data with the columns:
- observation #
- block code
- EE_seq denoting the code of the sequence of normoxia anoxia fluctuations during experimental evolution (see Figure S2)
- grandmother hatching environment, N for normoxia, A for anoxia
- mother hatching environment,  N for ormoxia, A for anoxia
- population internal lab code
- EEV_acronym, EEV wormbase lab code for each population
- EE, experimental regime, ANC for ancestral, REL for reliable, UNREL for unreliable
- cue, presence or absence of blue-light pulses in maternal generation during the assay (see Figure S1)
- growth_rate, absolute offspring fitness in anoxia (R_0)

cue_data.csv is superceded by cue_data
