# LCD-Composer
LCD-Composer is a ***l***ow-***c***omplexity ***d***omain ***compo***sition ***s***cann***er*** designed to identify simple or multifaceted low-complexity domains (LCDs) in protein sequences, described in ***(Cascarina et al. Genome Biology 2020) add link to reference when published***.

LCD-Composer identifies LCDs by calculating the amino acid composition and linear dispersion of amino acids at each position in a protein sequence using a scanning window of defined size. For a full description of how the algorithm works (complete with graphical representations of algorithm workflow and extensive parameter testing), see the publication cited above.
<br/>
## Basic Usage
$ python LCD-Composer.py Sequences_File Results_File [-a <amino_acids>] [-c <composition>] [-w window_size] [-d dispersion] [-i ignore_dispersion_threshold] [-v]

Identify LCDs based on amino acid composition and linear amino acid dispersion.

positional arguments:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Sequences_File &nbsp;&nbsp;&nbsp;&nbsp;The name of the file containing your protein sequences, in FASTA format.<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Results_File &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The name of the results file that you want to create and store the resulting LCD data.

required arguments:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-a <amino_acids> &nbsp;&nbsp;&nbsp;&nbsp;Amino acid(s) of interest (single letter abbreviation).<br/>

optional arguments:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-c \<composition> &nbsp;&nbsp;&nbsp;&nbsp;Percent composition threshold for amino acid(s) of interest (0-100).<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-w <window_size> &nbsp;&nbsp;&nbsp;Sliding window size.<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-d \<dispersion> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Linear dispersion threshold for amino acid(s) of interest (0.0-1.0).<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-i <ignore_dispersion_threshold> &nbsp;&nbsp;&nbsp;Threshold at which to ignore the linear dispersion parameter (0.0-1.0).<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-v \<verbose> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Output verbose data for each protein (per-position composition and dispersion values).<br/>

### Run LCD-Composer with Default Parameters
LCD-Composer is a Python3 script designed to be run as a stand-alone command-line application. To run LCD-Composer on your sequences of interest, download the LCD-Composer.py script and save to a location containing a FASTA file with your protein sequences of interest. Navigate to that location via command line, and run LCD-Composer with the following command (will use default parameters):

'''
python LCD-Composer.py Sequences_File Results_File -a <amino_acids>
'''

The "-a" flag is required, and allows you to specify the amino acid(s) that you are interested in for your domain search (see **Customizable Parameters** section below for complete description.

NOTE: Make sure to include the file extension in the command above for your file containing FASTA-formatted sequences. FASTA files will often have the file extension ".fa", ".fsa", or ".FASTA", but are sometimes also provided as plain-text files (.txt), which should still work with LCD-Composer. LCD-Composer is designed to output your results in a **t**ab-**s**eparated **v**alues (.tsv) file. This file type was chosen for two main reasons: 1) .tsv files can be easily parsed in downstream computational processing and avoids using comma-delimiters which are often present in FASTA headers, and 2) .tsv files can be opened by Microsoft Excel for the typical user. However, if Microsoft Excel is not set as the default program to open .tsv files, the file must be opened from *within* Excel (i.e. first open Excel, then open the results file from within Excel). Alternatively, you can first change your system settings to open .tsv files with Excel by default. Please note that Excel notoriously re-formats some gene/protein names to date formats automatically: however, this only occurs if the gene/protein name constitutes the entire FASTA record header in the sequence file, and only for a small number of genes/proteins (e.g. SEPT7, MARCH1, etc.).

## Detailed Usage and Customizable Parameters
Mulitple LCD-Composer parameters are customizable at runtime for more targeted CompoSer searches. These include:
1. Amino acid(s) of interest
2. Composition threshold(s)
3. Window size
4. Linear dispersion threshold
5. "Verbose" output

In the following sections, we illustrate the usage of each parameter.
<br/>
### Amino Acid(s) of Interest (-a)
LCD-Composer requires users to specify an amino acid or group(s) of amino acids of interest for each run. Amino acid(s) are specified using the "-a" flag followed by a space, followed by the single-letter abbreviation for the amino acid of interest. For example, a search for Q-rich regions (with other LCD-Composer parameters set to default), the command would be:

> \>python LCD-Composer.py Sequences_File Results_File -a Q

LCD-Composer also permits searches for domains enriched in a specific set of residues. For example, users can search for domains enriched in negatively charged residues (D and/or E) with the following command:

> \>python LCD-Composer.py Sequences_File Results_File -a DE

### Composition thresholds (-c)
The default composition threshold for LCD-Composer is 40% composition. This means that the specified amino acid(s) of interest must be at least 40% of a given sequence window for further consideration as a domain of interest. Alternative composition thresholds can be specified using the "-c" flag, followed by a space, followed by any value from 0-100 (inclusive). For example, to search for N-rich regions that are at least 60% N, the command would be:

> \>python LCD-Composer.py Sequences_File Results_File -a N -c 60

LCD-Composer also allows for "multifaceted" composition criteria, where distinct amino acids or groups of amino acids are assigned differnt thresholds. For example, users may be interested in domains that are at least 40% S ***and*** at least 20% A. This can be done by separating the amino acids with an underscore, while also separating distinct composition thresholds by an underscore. The command would be:

> \>python LCD-Composer.py Sequences_File Results_File -a S_A -c 40_20

This can also be performed with groups of amino acids, by separating groups of amino acids with the underscore delimiter. For example, a search for domains that are at least 50% Q and/or N ***and*** at least 15% Y, the command would be:

> \>python LCD-Composer.py Sequences_File Results_File -a QN_Y -c 50_15

Notice that any amino acids that are not separated by the underscore will be considered a group, and their combined composition for each window will be considered in the calculation.

### Window size (-w)
By default, LCD-Composer uses a sliding window size of 20 amino acids. To use an alternative window size, use the "-w" flag, followed by a space, followed by any positive integer value. For example, a search for long S-rich domains at least 60 residues in length would be:

> \>python LCD-Composer.py Sequences_File Results_File -a S -w 60

### Linear dispersion threshold (-d)
The linear dispersion parameter is a normalized measure of the dispersion of the amino acid(s) of interest within each window. Linear dispersion will always be a decimal value from 0-1, with higher values indicating greater linear dispersion (i.e. approaching perfect spacing of the amino acid(s) of interest within a given window sequence). By default, LCD-Composer uses a default linear dispersion threshold of 0.5, with the added caveat that the linear dispersion parameter is ignored if the composition of a window sequences exceeds the midpoint between the composition threshold and 100%. For example, at a composition threshold of 40%, the linear dispersion parameter will be ignored for windows with at least 70% composition corresponding to the amino acid(s) of interest.

To define a new linear dispersion threshold, e.g. 0.7, in a search for N-rich domains, the command would be:

> \>python LCD-Composer.py Sequences_File Results_File -a N -d 0.7

### "Verbose" output (-v)
Some users may be interested in the composition and linear dispersion values assigned to each position of the protein. This can be achieved by using the "-v" flag (with no other trailing arguments or characters). For example, the per-position composition and dispersion of Q residues can be assessed using the following command:

> \>python LCD-Composer.py Sequences_File Results_File -a Q -v

NOTE: verbose mode forces CompoSer to perform the complete set of calculations for each position in a protein. Consequently, verbose output runs will be slower than the default LCD-Composer, which only performs full calculations for identified LCDs.

### Combining customizable commands
All LCD-Composer commands can be used in combination for highly selective searches. Below are examples of combined commands, and a brief description of what each search is designed to accomplish:

> \>python LCD-Composer.py Sequences_File Results_File -a P -c 65 -w 50 -d 0.6
<br/>__(search for P-rich domains that are at least 65% P, at least 50aa long, and with moderate minimum dispersion of 0.6)__

> \>python LCD-Composer.py Sequences_File Results_File -a DE_KR -c 40_40 -w 30 -d 0.4
<br/>__(search for domains that are at least 40% D or E **and** at least 40% K or R (mixed charge domains), that are at least 30aa long, and with moderate minimum dispersion of 0.4)__

> \>python LCD-Composer.py Sequences_File Results_File -a G_RY -c 30_15 -w 60 -d 0.7
<br/>__(search for domains that are at least 30% G **and** at least 15% R or Y (often associated with RGG domains), that are at least 60aa long, and with relatively high minimum dispersion of 0.7)__

> \>python LCD-Composer.py Sequences_File Results_File -a FWY -c 25 -w 35 -d 0.8 -v
<br/>__(search for aromatic-rich domains that are at least 25% F, W, or Y; at least 35aa long; with high minimum dispersion of 0.8; and output the per-position values for each protein in "verbose" mode)__

## License info
LCD-Composer is subject to the terms of the GPLv3 license.
