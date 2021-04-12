# Reproducing Analyses Pertaining to Development of the Linear Dispersion Parameter

### Instructions
1. Download all files in the LinearDispersionDevelopment directory.
2. Navigate to appropriate folder via command line.
3. Run the following commands in-sequence to generate Fig S1 and to validate determination of the minimum and maximum linear dispersion for benchmark sequences from 5aa to 30aa in length (NOTE: each run must be completed before issuing the next command):

```
python LCD-Composer_Analyze_BenchmarkSequences.py 20aaBenchmarkSequences_LCD-Composer_RESULTS -a G -c 0 -d 0 -w 20 -l 20 -r 19
```

```
python plot_LinearDispersionDistributions_ViolinPlots.py
```

```
python make_LinearDispersionValidation_BatchFile.py
```

```
.\RUN_LinearDispersionValidation.bat
```

Run the following command to generate Fig S2:

```
python plot_DiscriminatoryEffect_of_LinearDispersion.py
```

Run the following commands in-sequence to generate Fig S3:

```
python build_Clustered_BenchmarkSequences.py
```

```
python make_ClusteredBenchmarkSequences_BatchFile.py
```

```
.\RUN_LCD-Composer_ClusteredBenchmarkSequences.bat
```

```
python plot_Breakpoints_for_ClusteredBenchmarkSequences.py
```
