# Reproducing Analyses Related to 7 Model Eukaryotic Organisms

### Instructions
1. Download LCD-Composer.py, as well as all files in the ModelEukaryoticOrganisms directory and place in the same folder.
2. Extract files from compressed folders in the same location as LCD-Composer.py
4. Navigate to appropriate folder via command line.
5. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

```
python write_LCD-Composer_BatchFile_ModelEukaryoticOrganisms_Analyses.py
```
```
.\Run_LCD-Composer_on_ModelEukaryotes_AAvaried.bat
```
```
python gather_ModelEukaryoticOrganisms_All_LCDs.py
```
```
python gather_LCD-ContainingProteins.py
```
```
python Run_GOanalyses_ModelEukaryoticOrganisms.py
```
```
python gather_ModelEukaryoticOrganisms_AllStatisticallySignificantGOterms.py
```
```
python Calculate_Homotypic_vs_Heterotypic_LCD_Overlap.py
```
```
python plot_ModelEukaryoticOrganisms_CrossOrganismGOtermFrequencies.py
```
```
python plot_ModelEukaryoticOrganisms_ProportionGOterms_STACKED_BARCHART.py
```

This series of commands generates Fig 6, Fig S7, Fig S12, and Tables S3-S6.
</br></br></br>
Run the following commands in-sequence to generate Fig S8:
```
python ModelOrganisms_LENGTH-WEIGHTED_ProteinSampling_GOtermAnalyses.py
```
```
python calculate_ModelEukaryoticOrganisms_LENGTH-WEIGHTED_ProteinSampling_Cross-Organism_GOfrequencies.py
```
```
python plot_ModelEukaryoticOrganisms_LENGTH-WEIGHTED_ProteinSampling_GOresults.py
```
</br>
GO term analysis involving all protein isoforms for each organism requires running LCD-Composer on a different version of each proteome. Run the following commands in-sequence to generate Fig S9:

```
python write_LCD-Composer_BatchFile_ModelEukaryoticOrganisms_Analyses_ALL-ISOFORMS.py
```
```
.\Run_LCD-Composer_on_ModelEukaryotes_AAvaried_ALL-ISOFORMS.bat
```
```
python gather_LCD-ContainingProteins_ALL-ISOFORMS.py
```
```
python Run_GOanalyses_ModelEukaryoticOrganisms_ALL-ISOFORMS.py
```
```
python plot_ModelEukaryoticOrganisms_ALL-ISOFORMS_CrossOrganismGOtermFrequencies.py
```
</br>
To run GO term analyses with homology-based GO evidence codes excluded, we first need to create gene association files (.gaf) with these homology codes excluded. Run the following commands in-sequence to generate Fig S10:

```
python ExcludedAnnotations_GAF_files.py
```
```
python Run_GOanalyses_ModelEukaryoticOrganisms_EXCLUDED-ANNOTS.py
```
```
python plot_ModelEukaryoticOrganisms_EXCLUDED-ANNOTS_CrossOrganismGOtermFrequencies.py
```
</br>
Run the following commands in-sequence to generate Fig S11 and Table S7 (NOTE: since the final script implements random sampling, actual data values in the resulting table may differ slightly from those in Table S9, but should be statistically comparable):

```
python get_proteome_GOterm_lists.py
```
```
python calculate_NumberAssociatedProts_from_GAF_files.py
```
```
python determine_OrgansimSpecific_and_LCDspecific_pval_cutoffs_for_GOtermSampling.py
```
```
python CrossSpecies_GOterm_Sampling_Estimation.py
```
</br>
Run the following commands in-sequence to generate Fig S6B and Fig S6C:

```
python LCD-Composer_Modified_VaryAA_and_Output_ComputationTimes.py
```
```
python plot_ModelOrganisms_Computation_Times.py
```
</br>
Run the following commands in-sequence to generate Table S9:

```
python make_ModelEukaryoticOrganisms_PrimaryAA-SecondaryAA_LCD-Composer_BatchFile.py
```
```
.\Run_LCD-Composer_ModelEukaryoticOrganisms_PrimaryAA-SecondaryAA.bat
```
```
python gather_LCDproteins_PrimaryAA-SecondaryAA_ModelEukaryoticOrganisms.py
```
```
python Run_GOterm_Analyses_PrimaryAA-SecondaryAA_ModelEukaryoticOrganisms.py
```
```
python Compare_PrimaryGOenrichment_vs_LCDsubclassingGOenrichment.py
```
</br>
Run the following commands in-sequence to generate Table S11:

```
python gather_MultiLCD_Proteins_ModelEukaryoticOrganisms.py
```
```
python Run_GOterm_Analyses_MultiLCD_Proteins_ModelEukaryoticOrganisms.py
```
```
python Compare_Primary_vs_Secondary_MultiLCD_GOenrichment.py
```
