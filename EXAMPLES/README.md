# Instructions

A set of test sequences has been prepared to demonstrate LCD-Composer usage and expected output.

1. Download the Test_Sequences.fasta file in the EXAMPLES folder.
2. Download LCD-Composer.py and place in the same folder as Test_Sequences.fasta.<br/>
3. Run any or all of the following commands:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;***-Basic search for Q-rich LCDs:***<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;python LCD-composer.py TestSequences.fasta TestSequences_Q_RESULTS -a Q<br/><br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;***-Search for extremely D-rich LCDs (at least 60% D):***<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;python LCD-composer.py TestSequences.fasta TestSequences_D_RESULTS -a D -c 60<br/><br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;***-Search for moderately E/K-rich LCDs (at least 50% E/K):***<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;python LCD-composer.py TestSequences.fasta TestSequences_EK_RESULTS -a EK -c 50<br/><br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;***-Search for long Q/N-rich LCDs with at least 10% Y and moderate dispersion (>0.6):***<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;python LCD-composer.py TestSequences.fasta TestSequences_QNY_RESULTS -a QN_Y -c 40_10 -w 60 -d 0.6<br/><br/>
4. Compare the resulting output files with the pre-computed results files stored in the EXAMPLES folder.
