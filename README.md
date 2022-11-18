## Spatially localized fMRI metrics as predictive and highly distinct state-independent fingerprints

This repository contains code to generate fingerprints and perform identification.
The details are described in the paper "Spatially localized fMRI metrics as predictive and highly distinct state-independent fingerprints".
The [preprint is available on biorxiv](https://www.biorxiv.org/content/10.1101/2021.08.03.454862v3).

## Code

Please run the `setup.sh` file which will get the DataLad dataset and spm12.
Note that you need [datalad](https://www.datalad.org/) available in your path.

The main analysis code is available in the `code` folder.
We also provide two convenience files `calculate_fingerprints.m` and `calculate_identification.m`.
The `calculate_fingerprints.m` will calculate fingerprints using the [fMRIPrep derivates of ds000115 dataset](https://github.com/OpenNeuroDerivatives/ds000115-fmriprep).
The script will download the DataLad dataset and use corresponding files to calculate several fingerprints as described
in the paper.
For demonstration only 10 subjects are used but this can be changed in the code.
The output will be stored in the `output` folder.
The `calculate_identification.m` script can be then used to calculate identification accuracies from those fingerprints.
For convenience we also provide mask and atlas files in the correct space.

This code has only been tested on a Linux machine.
