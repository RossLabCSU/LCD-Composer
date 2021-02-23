# Reproducing Analyses of Viruses Using LCD-Composer

### Instructions
1. Download LCD-Composer.py, as well as all files in the Viruses_Analyses directory
2. Download all viral proteomes from https://figshare.com/articles/dataset/Viruses_ProteinSequences/12942362
3. Extract proteome files from compressed file in the same location as LCD-Composer.py
4. Navigate to appropriate folder via command line.
5. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

>\>python make_Viruses_LCD-Composer_BatchFile.py

>\>RUN_LCD-Composer_Viruses_Batch.bat

>\>python calculate_Viruses_ProteomeCharacteristics_MeanMedianLCDcontent.py

>\>python plot_Viruses_Top_LCDcoverage_Proteomes.py

This series of commands generates data appearing in Fig 4, Fig 5D, and Table S1.
