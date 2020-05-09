# CompoSer
CompoSer is a ***compo***sition ***s***cann***er*** designed to identify simple or complex low-complexity domains (LCDs) in protein sequences, described in ***(Cascarina et al. eLife 2020) add link to reference when published***.

CompoSer identifies LCDs by calculating the amino acid composition and linear dispersion of amino acids at each position in a protein sequence using a scanning window of defined size. For a complete description of how the algorithm works (complete with graphical representations of algorithm workflow and extensive parameter testing), see the publication cited above.
<br/>
## Basic Usage
### Run CompoSer with Default Parameters
CompoSer is a Python3 script designed to be run as a stand-alone command-line application. To run CompoSer on your sequences of interest, download the CompoSer.py script and save to a location containing a FASTA file with your protein sequences of interest. Navigate to that location via command line, and run CompoSer with the following command (will use default parameters):

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file

NOTE: Make sure to include the file extension in the command above for your file containing FASTA-formatted sequences. FASTA files will often have the file extension ".fsa" or ".FASTA", but are sometimes also provided as plain-text files (.txt), which should still work with CompoSer. As specified above, CompoSer is designed to output your results in a **t**ab-**s**eparated **v**alues (.tsv) file. This file type was chosen for two main reasons: 1) .tsv files can be easily parsed in downstream computational processing and avoids using comma-delimiters which are often present in FASTA headers, and 2) .tsv files can be opened by Microsoft Excel for the typical user. However, if Microsoft Excel is not set as the default program to open .tsv files, the file must be opened from *within* Excel (i.e. first open Excel, then open the results file from within Excel). Alternatively, you can first change your system settings to open .tsv files with Excel by default.
<br/>
## Customizable Parameters
Mulitple CompoSer parameters are customizable at runtime for more targeted CompoSer searches. These include:
1. Amino acid(s) of interest
2. Composition threshold(s)
3. Window size
4. Linear dispersion threshold
5. "Verbose" output

In the following sections, we illustrate the usage of each parameter.
<br/>
### Amino Acid(s) of Interest (-a)
CompoSer requires users to specify an amino acid or group(s) of amino acids of interest for each run. Amino acid(s) are specified using the "-a" flag followed by a space, followed by the single-letter abbreviation for the amino acid of interest. For example, a search for Q-rich regions (with other CompoSer parameters set to default), the command would be:

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a Q

CompoSer also permits searches for domains enriched in a specific set of residues. For example, users can search for domains enriched in negatively charged residues (D and/or E) with the following command:

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a DE
<br/>
### Composition thresholds (-c)
The default composition threshold for CompoSer is 40% composition. This means that the specified amino acid(s) of interest must be at least 40% of a given sequence window for further consideration as a domain of interest. Alternative composition thresholds can be specified using the "-c" flag, followed by a space, followed by any value from 0-100 (inclusive). For example, to search for N-rich regions that are at least 60% N, the command would be:

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a N -c 60

CompoSer also allows for "complex" composition criteria, where distinct amino acids or groups of amino acids are assigned differnt thresholds. For example, users may be interested in domains that are at least 40% S ***and*** at least 20% A. This can be done by separating the amino acids with an underscore, while also separating distinct composition thresholds by an underscore. The command would be:

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a S_A -c 40_20

This can also be performed with groups of amino acids, by separating groups of amino acids with the underscore delimiter. For example, a search for domains that are at least 50% Q and/or N ***and*** at least 15% Y, the command would be:

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a QN_Y -c 50_15

Notice that any amino acids that are not separated by the underscore will be considered a group, and their combined composition for each window will be considered in the calculation.
<br/>
### Window size (-w)
By default, CompoSer uses a sliding window size of 20 amino acids. To use an alternative window size, use the "-w" flag, followed by a space, followed by any positive integer value. For example, a search for long S-rich domains at least 60 residues in length would be:

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a S -w 60
<br/>
### Linear dispersion threshold (-d)
The linear dispersion parameter is a normalized measure of the dispersion of the amino acid(s) of interest within each window. Linear dispersion will always be a decimal value from 0-1, with higher values indicating greater linear dispersion (i.e. approaching perfect spacing of the amino acid(s) of interest within a given window sequence). By default, CompoSer uses a default linear dispersion threshold of 0.5, with the added caveat that the linear dispersion parameter is ignored if the composition of a window sequences exceeds the midpoint between the composition threshold and 100%. For example, at a composition threshold of 40%, the linear dispersion parameter will be ignored for windows with at least 70% composition corresponding to the amino acid(s) of interest.

To define a new linear dispersion threshold, e.g. 0.7, in a search for N-rich domains, the command would be:

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a N -d 0.7
<br/>
### "Verbose" output (-v)
Some users may be interested in the composition and linear dispersion values assigned to each position of the protein. This can be achieved by using the "-v" flag (with no other trailing arguments or characters). For example, the per-position composition and dispersion of Q residues can be assessed using the following command:

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a Q -v

NOTE: verbose mode forces CompoSer to perform the complete set of calculations for each position in a protein. Consequently, verbose output runs will be slower than the default CompoSer, which only performs full calculations for identified LCDs.
<br/>
### Combining customizable commands
All CompoSer commands can be used in combination for highly selective searches. Below are examples of combined commands, and a brief description of what each search is designed to accomplish:

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a P -c 65 -w 50 -d 0.6
__(search for P-rich domains that are at least 65% P, at least 50aa long, and with moderate minimum dispersion of 0.6)__

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a DE_KR -c 40_40 -w 30 -d 0.4
__(search for domains that are at least 40% D or E **and** at least 40% K or R (mixed charge domains), that are at least 30aa long, and with moderate minimum dispersion of 0.4)__

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a G_RY -c 30_15 -w 60 -d 0.7
__(search for domains that are at least 30% G **and** at least 15% R or Y (often associated with RGG domains), that are at least 60aa long, and with relatively high minimum dispersion of 0.7)__

> \>python CompoSer.py -o your_results_file.tsv your_FASTA_sequences_file -a FWY -c 25 -w 35 -d 0.8 -v
__(search for aromatic-rich domains that are at least 25% F, W, or Y; at least 35aa long; with high minimum dispersion of 0.8; and output the per-position values for each protein in "verbose" mode)__

## License info
CompoSer is subject to the terms of the GPLv3 license.
