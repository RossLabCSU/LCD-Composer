# Reproducing Yeast-Specific LCD-Composer Analyses

### Instructions
1. Download LCD-Composer.py, as well as all files in the Yeast directory and place in the same folder.
2. Extract files from compressed folders in the same location as LCD-Composer.py
3. Navigate to appropriate folder via command line.
4. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):
</br>

```
python LCD-Composer_Modified_VaryWindowSize.py UP000002311_Scerevisiae_NoIsoforms.fasta
```
```
python LCD-Composer_Modified_VaryLinearDispersion.py UP000002311_Scerevisiae_NoIsoforms.fasta
```
```
python LCD-Composer_Modified_VaryCompositionThreshold.py UP000002311_Scerevisiae_NoIsoforms.fasta
```
```
python plot_Yeast_LCD_Lengths_ParametersVaried.py
```
```
python plot_Yeast_LCD-Composer_ComputationTimes.py
```
```
python plot_Yeast_TotalNumLCDs_Heatmaps_ParametersVaried.py
```
This series of commands generates all panels in Fig S4, all panels in Fig S5, and Fig S6A. NOTE: Computation times will, of course, vary based on computer hardware and usage.
</br></br></br>
Run the following commands in-sequence to generate Fig 8 and Table S8:

```
python Calculate_Background_X-rich_Window_Frequencies.py
```
```
python Plot_SecondaryAA_LCD_Enrichment.py
```
```
python Calculate_SecondaryAA_LCD_Enrichment_Statistics.py
```
```
python Apply_Bonferroni_Correction_SecondaryAA_LCD_Enrichment_Statistics.py
```
</br>
Initial GO term analysis of the 20 main LCD classes is performed by running the first 4 commands indicated under step #4 in the instructions in the ModelEukaryoticOrganisms directory. Once these commands have been completed, copy the 19 resulting "Scerevisiae_X_GO_RESULTS.tsv" files and the "Scerevisiae_X_LCD-containing_proteins.txt" files (X represents each of the amino acids...the file for W is never generated because there are no W-rich LCDs in the yeast proteome by these search criteria) into the same folder containing the downloaded Yeast files. Then run the following commands in-sequence:

```
python make_Yeast_PrimaryAA-SecondaryAA_LCD-Composer_BatchFile.py
```
```
.\Run_LCD-Composer_Yeast_PrimaryAA-SecondaryAA.bat
```
```
python gather_LCDproteins_PrimaryAA-SecondaryAA.py
```
```
python Run_GOterm_Analyses_PrimaryAA-SecondaryAA_Yeast.py
```
```
python Compare_PrimaryGOenrichment_vs_LCDsubclassificationGOenrichment.py
```
```
python plot_LCDsubclassification_GOresults_stats.py
```
```
python Compare_lnORs_Primary_vs_Secondary_LCDs.py
```
This series of commands generates Fig 9C and Table S10.
</br></br></br>
Removal of proteins with high homology from LCD sets prior to analysis (pertaining to Fig S15) requires multiple sequence alignment using the Clustal Omega server submitted via an API script. Run the following commands in-sequence to generate Fig S15A:

```
python gather_PrimaryLCD_Proteins.py
```
```
python make_PrimaryLCD_ProteinSequences_FastaFiles.py
```
```
python make_ClustalW_BatchFile_PrimaryLCDs.py
```
```
.\Run_Yeast_PrimaryLCDs_Alignments_ClustalW_BatchFile.bat
```
```
python plot_ClustalW_PrimaryLCDs_PercIdentity_Boxplots.py
```
</br>
Run the following commands in-sequence to generate Fig S15B:

```
python make_LCDsubclassification_ProteinSequences_FastaFiles.py
```
```
python make_ClustalW_BatchFile_LCDsubclassification.py
```
```
.\Run_Yeast_LCDsubclassification_Alignments_ClustalW_BatchFile.bat
```
```
python get_HighHomology_Prots_LCDsubclassification.py
```
```
python remove_HighHomologyProts_SecondaryAA_Prots.py
```
```
python Run_LCDsubclassification_GOtermAnalyses_HighHomologyProtsRemoved.py
```
```
python Compare_Primary_vs_Secondary_GOenrichment_80percHomologyProtsRemoved.py
```
```
python plot_Primary_vs_Secondary_GOenrichment_80percHomologyProtsRemoved.py
```
</br>
Run the following commands in-sequence to generate Fig 10B, Fig 10C, and Table S12:

```
python gather_MultiLCD_Proteins_and_plot_MultiLCD_Frequencies_Yeast.py
```
```
python Run_GOanalyses_Secondary_AA_MultiLCD_proteins.py
```
```
python Compare_Primary_vs_Secondary_GOenrichment_MultiLCD_Proteins.py
```
```
python plot_LCDsubclassification_GOresults_stats_MultiLCD_Proteins.py
```
```
python Compare_lnORs_MultiLCD_Proteins.py
```
</br>
Run the following commands in-sequence to generate Fig S15C:

```
python make_MultiLCD_ProteinSequences_FastaFiles.py
```
```
python make_ClustalW_BatchFile_MultiLCD_Proteins.py
```
```
.\Run_Yeast_MultiLCD_Alignments_ClustalW_BatchFile.bat
```
```
python get_HighHomology_Prots_MultiLCD.py
```
```
python remove_HighHomologyProts_MultiLCD_Prots.py
```
```
python Run_MultiLCD_GOtermAnalyses_HighHomologyProtsRemoved.py
```
```
python Compare_MultiLCD_GOenrichment_80percHomologyProtsRemoved.py
```
```
python plot_MultiLCD_GOenrichment_80percHomologyProtsRemoved.py
```
</br>
The following command generates all heatmaps appearing in Fig 7, Fig S13, and Fig S14:

```
python plot_LCD_Composition_Clustermaps_Yeast.py
```
