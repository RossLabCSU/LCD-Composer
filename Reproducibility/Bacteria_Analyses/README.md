# Reproducing Analyses of Bacteria Using LCD-Composer

### Instructions
1. Download LCD-Composer.py, as well as all files in the Bacteria_Analyses directory
2. Download all bacterial proteomes from https://figshare.com/articles/dataset/Bacteria_ProteinSequences_Part1/12939812 and https://figshare.com/articles/dataset/Bacteria_ProteinSequences_Part2/12939980
3. Extract proteome files from compressed file in the same location as LCD-Composer.py
4. Download the "TaxonID_to_Organism.dat" file from the main "Reproducibility" directory and place in the same location as the proteome files and LCD-Composer.py
5. Navigate to appropriate folder via command line.
6. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

>\>python make_Bacteria_LCD-Composer_BatchFile.py

>\>RUN_LCD-Composer_Bacteria_Batch.bat

>\>python calculate_Bacteria_ProteomeCharacteristics_MeanMedianLCDcontent.py

>\>python plot_Bacteria_Top_LCDcoverage_Proteomes.py

This series of commands generates data appearing in Fig 4, Fig 5B, and Table S1.
