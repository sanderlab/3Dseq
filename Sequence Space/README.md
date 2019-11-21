
# Sequence Space Scripts

Requires `sequence_space_data.zip` from https://evcouplings.org/3Dseq

1. `Load Taxonomy Labels from Uniprot.ipynb`
   This notebook contains code used to download taxonomy labels sequences. It was
   used to determine the class labels and color individual sequences in the
   sequence space figure.
     - Input data = fasta file
     - Output data = a tab separated file of sequence names and their taxonomy labels
2. `executevae_PSE1.py and execute_AAC6.py`
   These scripts generate a non-linear latent variable model given a subsampling (5000)
   of natural homologs of PSE1 or AAC6. The model is used to project down from
   N-dimensional space into two dimensional space - see the jupyter notebook
   `Generate sequence landscape figure.ipynb`
    - Input data = an alignment containing natural homologs of PSE1 or AAC6
    - Output data = a non-linear latent variable model
    - Note - this model framework is based on work described in Riesselman et al., 
             Nature Methods, 2018
3. `train.py, model.py, helper.py`
   These scripts are taken from DeepSequence and are used to train the non-linear
   latent variable model. (see https://github.com/debbiemarkslab/DeepSequence)
4. `Generate sequence landscape figure.ipynb`
   This notebook contains the code that generates the sequence space figure.

