# Reproducing Analyses Pertaining to Development of the Linear Dispersion Parameter

### Instructions
1. Download all files in the LinearDispersionDevelopment directory.
2. Navigate to appropriate folder via command line.
3. Run the following commands in-sequence to generate Fig S1 (NOTE: each run must be completed before issuing the next command):

>\>python generate_20aaBenchmarkSequences.py

>\>python python LCD-Composer_Modified_to_Omit_LCD-Trimming.py 20aaBenchmarkSequences.FASTA 20aaBenchmarkSequences_LCD-Composer_RESULTS -a G -c 0 -d 0

>\>python plot_LinearDispersionDistributions_ViolinPlots.py

Run the following command to generate Fig S2:

>\>python plot_DiscriminatoryEffect_of_LinearDispersion.py

Run the following commands in-sequence to generate Fig S3:

>\>python build_Clustered_BenchmarkSequences.py

>\>python make_ClusteredBenchmarkSequences_BatchFile.py

>\>RUN_LCD-Composer_ClusteredBenchmarkSequences

>\>python plot_Breakpoints_for_ClusteredBenchmarkSequences.py
