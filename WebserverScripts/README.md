# LCD-Composer Command-line Scripts
The scripts in this directory are command-line adaptations of the options available on the LCD-Composer webserver (http://lcd-composer.bmb.colostate.edu). Below is a brief description of usage and arguments available for each script. For more detailed description of how to use LCD-Composer command-line arguments, please refer to the README file in the main LCD-Composer directory or [this Help page on the LCD-Composer webserver](http://lcd-composer.bmb.colostate.edu:9000/help). You can also read our short paper describing the LCD-Composer webserver (ADD LINK WHEN PUBLISHED).

**Skip to:**<br />
[Option 1: Standard LCD-Composer search](#option-1-standard-lcd-composer-search)<br />
[Option 2: LCD similarity search](#option-2-lcd-similarity-search)<br />
[Option 3: Single-protein LCD search](#option-3-single-protein-lcd-search)<br />
[Option 4: LCD enrichment analysis](#option-4-lcd-enrichment-analysis)<br />
[Option 5: LCD co-occurrence analysis](#option-5-lcd-co-occurrence-analysis)<br />
[Automated GO-term analysis](#automated-go-term-analysis)

## Dependencies
| Script | Dependencies* | Required Files |
| --- | --- | --- |
| Option1_LCDcomposer.py | *None* &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;| Sequences_File (FASTA format) |
| Option2_LCDsimilaritySearch.py | *None* | Sequences file (FASTA format); Query LCD sequence file (FASTA format) |
| Option3_SingleProtein.py | Matplotlib v3.1.1<br/>Seaborn v0.9.0 | Sequence file (FASTA format, only 1 sequence) |
| Option4_LCD_EnrichmentTest.py | Matplotlib v3.1.1<br/>Seaborn v0.9.0<br/>Scipy v1.3.1<br/>mpmath v1.1.0 | Sequences file (FASTA format); UniProt IDs file (one UniProt ID per line) |
| Option5_LCD_CoOccurrence.py | *None* | Jobfile IDs (i.e. numerical IDs of results files) |
| Run_GOterm_analysis.py | GOATOOLS v1.0.2<br/>Scipy v1.3.1 | Jobfile ID (i.e. numerical ID of results file); Sequences file (FASTA format, from UniProt; defines "background" protein set); Gene annotation file (GAF) |

*\*Compatibility with other versions of dependencies is possible but not guaranteed. All scripts were tested and functional with Python v3.7.4.*
<br/><br/>
# Option 1: Standard LCD-Composer search
**Basic usage:**

    python Option1_LCDcomposer.py Sequences_File [-a AMINO_ACIDS] [-c COMPOSITION] [-w WINDOW_SIZE] [-d DISPERSION] [-i IGNORE_DISPERSION_THRESHOLD] [-m METHOD] [-if ISOFORMS_FILE] [-r]

***positional arguments:***<br/>
| Positional Argument | Description |
| --- | --- |
| Sequences_File | The name of the file containing the protein sequences you wish to search for LCDs (in FASTA format). |
<br />

***required arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -a AMINO_ACIDS | --amino_acids AMINO_ACIDS | Amino acid(s) of interest (single letter abbreviation, with a single underscore between each amino acid group when desired). This becomes optional if the -r flag is used. |
<br />

***optional arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -c COMPOSITION | --composition COMPOSITION | Percent composition threshold for amino acid(s) of interest (0-100). |
| -w WINDOW_SIZE | --window_size WINDOW_SIZE | Sliding window size. |
| -d DISPERSION | --dispersion DISPERSION | Linear dispersion threshold for amino acid(s) of interest (0.0-1.0). |
| -i IGN_DISP_THRESHOLD | --ignore_dispersion_threshold IGN_DISP_THRESHOLD | Threshold at which to ignore the linear dispersion parameter (0.0-1.0). |
| -m METHOD | --dispersion_method METHOD | Linear dispersion method (accepted arguments are "New" and "Original": see [this Help page](http://lcd-composer.bmb.colostate.edu:9000/help) for more details). |
| -if ISOFORMS_FILE | --isoforms_file ISOFORMS_FILE | Specify another FASTA-formatted file containing additional isoforms to include in the search if desired. |
| -r | --run_autodetect | Run "auto-detect" mode (i.e. a separate LCD search for each of the 20 canonical amino acids). This will override any amino acid specified. |
<br />

# Option 2: LCD similarity search
**Basic usage:**

    python Option2_LCDsimilaritySearch.py Sequences_File Query_Sequence_File [-n NUMBER_OF_FEATURES] [-i IGNORE_DISPERSION_THRESHOLD] [-if ISOFORMS_FILE]

***positional arguments (order matters!):***<br/>
| Positional Argument | Description |
| --- | --- |
| Sequences_File | The name of the file containing the protein sequences you want to search for LCDs (in FASTA format). |
| Query_Sequence_File | The name of the file containing the LCD sequence for which you want to find similar LCDs (in FASTA format. |
<br />

***optional arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -n NUMBER_OF_FEATURES | --n_features NUMBER_OF_FEATURES | The number of defining features (amino acids) of your query LCD sequence that you wish to use as search parameters to find similar LCDs. |
| -i IGN_DISP_THRESHOLD | --ignore_dispersion_threshold IGN_DISP_THRESHOLD | Threshold at which to ignore the linear dispersion parameter (0.0-1.0). |
| -if ISOFORMS_FILE | --isoforms_file ISOFORMS_FILE | Specify another FASTA-formatted file containing additional isoforms to include in the search if desired. |
<br />

# Option 3: Single-protein LCD search
**Basic usage:**

    python Option3_SingleProtein.py Sequences_File [-a AMINO_ACIDS] [-c COMPOSITION] [-w WINDOW_SIZE] [-d DISPERSION] [-i IGNORE_DISPERSION_THRESHOLD] [-m METHOD] [-if ISOFORMS_FILE] [-r] [-p] [-ir IMAGE_RESOLUTION] [-iw IMAGE_WIDTH] [-ih IMAGE_HEIGHT] [-if IMAGE_FILETYPE] [-tl THRESHOLD_LINE] [-cp COLOR_PALETTE] [-pa PLOTTING_AAS]
    
***positional arguments:***<br/>
| Positional Argument | Description |
| --- | --- |
| Sequences_File | The name of the file containing the protein sequence you wish to search for LCDs (in FASTA format). |
<br />

***required arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -a AMINO_ACIDS | --amino_acids AMINO_ACIDS | Amino acid(s) of interest (single letter abbreviation, with a single underscore between each amino acid group when desired). This becomes optional if the -r flag is used. |
<br />

***optional arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -c COMPOSITION | --composition COMPOSITION | Percent composition threshold for amino acid(s) of interest (0-100). |
| -w WINDOW_SIZE | --window_size WINDOW_SIZE | Sliding window size. |
| -d DISPERSION | --dispersion DISPERSION | Linear dispersion threshold for amino acid(s) of interest (0.0-1.0). |
| -i IGN_DISP_THRESHOLD | --ignore_dispersion_threshold IGN_DISP_THRESHOLD | Threshold at which to ignore the linear dispersion parameter (0.0-1.0). |
| -m METHOD | --dispersion_method METHOD | Linear dispersion method (accepted arguments are "New" and "Original": see [this Help page](http://lcd-composer.bmb.colostate.edu:9000/help) for more details). |
| -r | --run_autodetect | Run "auto-detect" mode (i.e. a separate LCD search for each of the 20 canonical amino acids). This will override any amino acid specified. |
| -p | --plot | Generate a composition plot for your protein. Remaining parameters apply only when the -p/--plot option is specified. |
| -ir IMAGE_RESOLUTION | --img_resolution IMAGE_RESOLUTION | Resolution of plot image. |
| -iw IMAGE_WIDTH | --img_width IMAGE_WIDTH | Width of plot image in mm. |
| -ih IMAGE_HEIGHT | --img_height IMAGE_HEIGHT | Height of plot image in mm. |
| -if IMAGE_FILETYPE | --img_filetype IMAGE_FILETYPE | Filetype (i.e. file extension) of plot. |
| -tl THRESHOLD_LINE | --threshold_line THRESHOLD_LINE | Composition threshold at which to draw a horizontal threshold line. |
| -cp COLOR_PALETTE | --color_palette COLOR_PALETTE | Color palette for plotting LCD categories. Must be a series of hex codes separated by underscores and *without* the customary hashtag as the first character. |
| -pa PLOTTING_AAS | --plotting_aas PLOTTING_AAS | Amino acid(s) for plotting (if different from amino acids used in the LCD search). Each AA group should be separated by an underscore. |
<br />

# Option 4: LCD enrichment analysis
**Basic usage:**

    python Option4_LCD_EnrichmentTest.py Sequences_File [-a AMINO_ACIDS] [-c COMPOSITION] [-w WINDOW_SIZE] [-d DISPERSION] [-i IGNORE_DISPERSION_THRESHOLD] [-m METHOD] [-if ISOFORMS_FILE] [-r] [-iv]

***positional arguments:***<br/>
| Positional Argument | Description |
| --- | --- |
| Sequences_File | The name of the file containing the protein sequences you wish to search for LCDs (in FASTA format). |
<br />

***required arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -a AMINO_ACIDS | --amino_acids AMINO_ACIDS | Amino acid(s) of interest (single letter abbreviation, with a single underscore between each amino acid group when desired). This becomes optional if the -r flag is used. |
<br />

***optional arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -c COMPOSITION | --composition COMPOSITION | Percent composition threshold for amino acid(s) of interest (0-100). |
| -w WINDOW_SIZE | --window_size WINDOW_SIZE | Sliding window size. |
| -d DISPERSION | --dispersion DISPERSION | Linear dispersion threshold for amino acid(s) of interest (0.0-1.0). |
| -i IGN_DISP_THRESHOLD | --ignore_dispersion_threshold IGN_DISP_THRESHOLD | Threshold at which to ignore the linear dispersion parameter (0.0-1.0). |
| -m METHOD | --dispersion_method METHOD | Linear dispersion method (accepted arguments are "New" and "Original": see [this Help page](http://lcd-composer.bmb.colostate.edu:9000/help) for more details). |
| -if ISOFORMS_FILE | --isoforms_file ISOFORMS_FILE | Specify another FASTA-formatted file containing additional isoforms to include in the search if desired. |
| -r | --run_autodetect | Run "auto-detect" mode (i.e. a separate LCD search for each of the 20 canonical amino acids). This will override any amino acid specified. |
| -iv | --impute_value | Impute biased frequency value(s) for enrichment test if no LCDs were found. |
<br />

# Option 5: LCD co-occurrence analysis
**Basic usage:**

    python Option5_LCD_CoOccurrence.py Jobfile_IDs

***positional arguments:***<br/>
| Positional Argument | Description |
| --- | --- |
| Jobfile_IDs | The randomly generated numeric IDs of the results files you wish to compare. These numeric IDs are at the start of the file name for each file generated by Options1-4, and are also listed as runtime parameters inside each results file. |
<br />

# Automated GO-term analysis
**Basic usage:**

    python Run_GOterm_Analysis.py Jobfile_ID Sequences_File GAF_File

***positional arguments (order matters!):***<br/>
| Positional Argument | Description |
| --- | --- |
| Jobfile_ID | The randomly generated numeric ID of the results file you wish to perform a GO-term analysis on. These numeric IDs are at the start of the file name for each file generated by Options1-5, and are also listed as runtime parameters inside each results file. |
| Sequences_File | The name of the file containing the protein sequences you wish to search for LCDs (in FASTA format). |
| GAF_File | The name of the gene annotation file (GAF) to use for the GO-term analysis. Appropriate GAF files are available in the "GAF_files" directory |
<br />

## License info
LCD-Composer and its derivatives in this directory are subject to the terms of the GPLv3 license.
