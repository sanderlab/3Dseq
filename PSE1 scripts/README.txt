For PSE1 experimental evolution, scripts for filtering and analyzing fastq sequences.

1) ProcessPacBioSequel_PSE-1.m
MATLAB script processes the fastq reads, filtering for flanking sequences, Q-score, full-length, and removes possible contaminating sequences
Inputs - PacBio sequencing results in fastq format (P20.fastq)
Output - fasta file with processes sequences 

2) scriptMaskContactMap.m (Frank Poelwijk, in "AAC6 scripts" folder) can be adapted for PSE1 for inferred contact filtering
Inputs - ContactListPSE1.mat; L=266 

