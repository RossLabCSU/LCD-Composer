# Reproducing Analyses of Archaea Using LCD-Composer

### Instructions
1. Download LCD-Composer.py, as well as all files in the Archaea_Analyses directory
2. Download all archaeal proteomes from https://figshare.com/articles/dataset/Archaea_ProteinSequences/12937637
3. Extract proteome files from compressed file in the same location as LCD-Composer.py
4. Navigate to appropriate folder via command line.
5. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

>\>python make_Archaea_LCD-Composer_BatchFile.py

>\>RUN_LCD-Composer_Archaea_Batch.bat

>\>python calculate_ProteomeCharacteristics_MeanMedianLCDcontent.py

>\>python plot_Top_LCDcoverage_Proteomes.py

This series of commands generates data appearing in Fig 4, Fig 5A, and Table S1.