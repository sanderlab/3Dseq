
# Structural Variability Scripts

Requires `structural_variability_data.zip` from https://evcouplings.org/3Dseq

1. `n_choose_2_rmsd.m`
   This script takes a set of models, executes a superimposition between all pairwise
   models, determines the minimum atom distance between each residue pair, and outputs
   a list of all pairwise distances for each residue.

2. `Structural Variability Visualization.ipynb`
   This script uses the nchoose2 output files (see above) to calculate RMSD for each
   position and output those RMSD values to a file. This file can than be used in
   pymol to highlight the structural variability using the putty visualization (sausage)
   derived by applying the loadBfacts script (doi 10.6084/m9.figshare.1176991.v1).
