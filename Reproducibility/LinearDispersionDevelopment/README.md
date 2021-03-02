# Reproducing Analyses Pertaining to Development of the Linear Dispersion Parameter

### Instructions
1. Download LCD-Composer.py, as well as all files in the LinearDispersionDevelopment directory and place in the same folder.
2. Navigate to appropriate folder via command line.
3. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

>\>python generate_20aaBenchmarkSequences.py

>\>python python LCD-Composer_Modified_to_Omit_LCD-Trimming.py 20aaBenchmarkSequences.FASTA 20aaBenchmarkSequences_LCD-Composer_RESULTS -a G -c 0 -d 0

>\>python .py
