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

This series of commands generates all panels in Fig S5.

Run the following commands in-sequence to generate Fig S8:
