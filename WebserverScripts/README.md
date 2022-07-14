# LCD-Composer Webserver Command-line Scripts
The scripts in this directory are command-line adaptations of the options available on the LCD-Composer webserver (http://lcd-composer.bmb.colostate.edu). Below is a brief description of usage and arguments available for each script. For more detailed description of how to use LCD-Composer command-line arguments, please refer to the README file in the main LCD-Composer directory.

**Skip to:**<br />
[Option 1: Standard LCD-Composer search](#option-1-standard-lcd-composer-search)<br />
[Option 2: LCD similarity search](#option-2-lcd-similarity-search)<br />
[Option 3: Single-protein LCD search](#option-3-single-protein-lcd-search)<br />
[Option 4: LCD enrichment analysis](#option-4-lcd-enrichment-analysis)<br />
[Option 5: LCD co-occurrence analysis](#option-5-lcd-co-occurrence-analysis)<br />
[Automated GO-term analysis](#-automated-go-term-analysis)

# Option 1: Standard LCD-Composer search
**Basic usage:**

    python Option1_LCDcomposer.py Sequences_File [-a AMINO_ACIDS] [-c COMPOSITION] [-w WINDOW_SIZE] [-d DISPERSION] [-i IGNORE_DISPERSION_THRESHOLD] [-m METHOD] [-if ISOFORMS_FILE] [-r]

***positional arguments:***<br/>
| Positional Argument | Description |
| --- | --- |
| Sequences_File | The name of the file containing the protein protein sequences you wish to search for LCDs (in FASTA format). |
<br />

***required arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -a AMINO_ACIDS | --amino_acids AMINO_ACIDS | Amino acid(s) of interest (single letter abbreviation, with a single underscore between each amino acid group when desired). |
<br />

***optional arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -c COMPOSITION | --composition COMPOSITION | Percent composition threshold for amino acid(s) of interest (0-100). |
| -w WINDOW_SIZE | --window_size WINDOW_SIZE | Sliding window size. |
| -d DISPERSION | --dispersion DISPERSION | Linear dispersion threshold for amino acid(s) of interest (0.0-1.0). |
| -i IGN_DISP_THRESHOLD | --ignore_dispersion_threshold IGN_DISP_THRESHOLD | Threshold at which to ignore the linear dispersion parameter (0.0-1.0). |
| -m METHOD | --dispersion_method METHOD | Linear dispersion method (accepted arguments are "New" and "Original": see http://lcd-composer.bmb.colostate.edu/help for more details). |
| -r | --run_autodetect | Run "auto-detect" mode (i.e. a separate LCD search for each of the 20 canonical amino acids). This will override any amino acid specified. |
| -if ISOFORMS_FILE | --isoforms_file ISOFORMS_FILE | Specify another FASTA-formatted file containing additional isoforms to include in the search if desired. |
<br />

# Option 2: LCD similarity search
**Basic usage:**

    python Option2_LCDsimilaritySearch.py Sequences_File Query_Sequence_File [-n NUMBER_OF_FEATURES] [-i IGNORE_DISPERSION_THRESHOLD] [-if ISOFORMS_FILE]

***positional arguments (order matters!):***<br/>
| Positional Argument | Description |
| --- | --- |
| Sequences_File | The name of the file containing the protein protein sequences you want to search for LCDs (in FASTA format). |
| Query_Sequence_File | The name of the file containing the LCD sequence for which you want to find similar LCDs (in FASTA format. |
<br />

***optional arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -n NUMBER_OF_FEATURES | --n_features NUMBER_OF_FEATURES | The number of defining features (amino acids) of your query LCD sequence that you wish to use as search parameters to find similar LCDs. |
| -i IGN_DISP_THRESHOLD | --ignore_dispersion_threshold IGN_DISP_THRESHOLD | Threshold at which to ignore the linear dispersion parameter (0.0-1.0). |
| -if ISOFORMS_FILE | --isoforms_file ISOFORMS_FILE | Specify another FASTA-formatted file containing additional isoforms to include in the search if desired. |
<br />

# Option 3: Single-protein lcd search

# Option 4: LCD enrichment analysis

# Option 5: LCD co-occurrence analysis

# Automated GO-term analysis

## License info
LCD-Composer and its derivatives in this directory are subject to the terms of the GPLv3 license.
