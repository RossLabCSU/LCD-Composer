# Reproducing Analyses Pertaining to Development of the Linear Dispersion Parameter

### Instructions
1. Download all files in the LinearDispersionDevelopment directory.
2. Navigate to appropriate folder via command line.
3. Run the following commands in-sequence to generate Fig S1 (NOTE: each run must be completed before issuing the next command):

>\>python generate_BenchmarkSequences.py

>\>python LCD-Composer_Modified_to_Omit_LCD-Trimming.py 20aaBenchmarkSequences.FASTA 20aaBenchmarkSequences_LCD-Composer_RESULTS -a G -c 0 -d 0

>\>python plot_LinearDispersionDistributions_ViolinPlots.py

Run the following command to generate Fig S2:

>\>python plot_DiscriminatoryEffect_of_LinearDispersion.py

Run the following commands in-sequence to generate Fig S3:

>\>python build_Clustered_BenchmarkSequences.py

>\>python make_ClusteredBenchmarkSequences_BatchFile.py

>\>RUN_LCD-Composer_ClusteredBenchmarkSequences

>\>python plot_Breakpoints_for_ClusteredBenchmarkSequences.py

If you wish to generate and analyze benchmark sequences of a different length than the standard 20aa, you can use optional arguments to specify the sequence length (-l) and mutation range (-r). Mutation range is the maximum number of "mutations", which represents the maximum number of X's in the otherwise polyG sequences. HOWEVER, note that file sizes become quickly extremely large as sequence length is increased. For example, the file size generated using a sequence length =30 and mutation range =15 will generate all necessary 30aa sequences to fully test the linear dispersion parameter, but the file generated is >35GB in size. This can be performed using the following command:

>\>python generate_BenchmarkSequences.py -l 30 -r 15

The mutation range (-r) =15 only generates sequences up to 50% composition for the two amino acids, which is the minimum necessary to test the linear dispersion parameter since reciprocal sequences will result in identical linear dispersion values. However, if you wish to generate the complete set of possible 30aa sequences use the "-r 29" flag, and the resulting file size should be >62GB. To generate the full set of 30aa sequences, run:

>\>python generate_BenchmarkSequences.py -l 30 -r 29

Then to analyze the resulting file using LCD-Composer, run:

>\>python LCD-Composer_Modified_to_Omit_LCD-Trimming.py 30aaBenchmarkSequences.FASTA 30aaBenchmarkSequences_LCD-Composer_RESULTS -a G -c 0 -d 0
