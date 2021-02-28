# Reproducing Yeast-Specific LCD-Composer Analyses

### Instructions
1. Download LCD-Composer.py, as well as all files in the Yeast directory and place in the same folder.
2. Extract files from compressed folders in the same location as LCD-Composer.py
4. Navigate to appropriate folder via command line.
5. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

>\>python LCD-Composer_Modified_VaryWindowSize.py UP000002311_Scerevisiae_NoIsoforms.fasta

>\>python LCD-Composer_Modified_VaryLinearDispersion.py UP000002311_Scerevisiae_NoIsoforms.fasta

>\>python LCD-Composer_Modified_CompositionThreshold.py UP000002311_Scerevisiae_NoIsoforms.fasta

>\>python plot_Yeast_LCD_Lengths_ParametersVaried.py

>\>python plot_Yeast_LCD-Composer_ComputationTimes.py

>\>python plot_Yeast_TotalNumLCDs_Heatmaps_ParametersVaried.py

This series of commands generates all panels in Fig S4, all panels in Fig S5, and Fig S6A. NOTE: Computation times will, of course, vary based on computer hardware and usage.

Initial GO term analysis of the 20 main LCD classes is performed by running the first 4 commands indicated under step #4 in the instructions in the ModelEukaryoticOrganisms directory. Once these commands have been completed, copy the 19 resulting "Scerevisiae_X_GO_RESULTS.tsv" files (X represents each of the amino acids...the file for W is never generated because there are no W-rich LCDs in the yeast proteome by these search criteria) into the same folder containing the downloaded Yeast files. Then run the following commands in-sequence:

>\>python make_Yeast_PrimaryAA-SecondaryAA_LCD-Composer_BatchFile.py

>\>Run_LCD-Composer_Yeast_PrimaryAA-SeoncdaryAA.bat

>\>python gather_LCDproteins_PrimaryAA-SecondaryAA.py

>\>python Run_GOterm_Analyses_PrimaryAA-SecondaryAA_Yeast.py

>\>python Compare_PrimaryGOenrichment_vs_LCDsubclassificationGOenrichment.py

>\>python plot_LCDsubclassification_GOresults_stats.py

>\>python Compare_lnORs_Primary_vs_Secondary_LCDs.py

This series of commands generates Fig 9C and Table S12.



