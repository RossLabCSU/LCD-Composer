# Reproducing Analyses Related to 7 Model Eukaryotic Organisms

### Instructions
1. Download LCD-Composer.py, as well as all files in the ModelEukaryoticOrganisms directory and place in the same folder.
2. Extract files from compressed folders in the same location as LCD-Composer.py
4. Navigate to appropriate folder via command line.
5. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

>\>python write_LCD-Composer_BatchFile_ModelEukaryoticOrganisms_Analyses.py

>\>Run_LCD-Composer_on_ModelEukaryotes_AAvaried.bat

>\>python gather_LCD-ContainingProteins.py

>\>python Run_GOanalyses_ModelEukaryoticOrganisms.py

>\>python Calculate_Homotypic_vs_Heterotypic_LCD_Overlap.py

>\>python ModelOrganisms_Non-Weighted_ProteinSampling_GOtermAnalyses.py

>\>python calculate_ModelEukaryoticOrganisms_Non-Weighted_ProteinSampling_Cross-Organism_GOfrequencies.py

>\>python plot_ModelEukaryoticOrganisms_Non-Weighted_ProteinSampling_GOresults.py

This series of commands generates Fig 6, Fig S7, Table S4, Table S6, Table S7, Table S8, and all data appearing in Table S3 and Table S5.
