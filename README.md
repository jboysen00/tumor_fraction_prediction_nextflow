# Tumor Fraction Prediction Nextflow
Joanne Boysen

July 8th, 2025

## Overview
This repository provides a Nextflow pipeline for tumor fraction prediction. Starting from BAM files aligned to the HG19 human genome reference, the pipeline generates BED files and tumor fraction estimates based on the HG38 reference. The workflow is fully containerized using Docker for reproducibility.

## Instructions for running the pipeline
### Build the Docker image
You must have Nextflow and Docker installed on your system to run this pipeline. Once Docker is installed, you can build the Docker image using the following script:
```
docker buildx build --platform linux/amd64 --no-cache -t tumor_fraction_prediction .
```

### Running the pipeline with example input data
If you wish to test the pipeline with example input data from `input_sample_data_200K.tar`, then you can run the pipeline with the following script:
```
nextflow run main.nf
```
If you would like an execution report with resource usage analysis to be generated when running the pipeline you can add `-with-report`:
```
nextflow run main.nf -with-report
```
An output table with each sample's tumor fraction can be found at `test_project/output/sample_tumor_fractions.tsv`

### Running the pipeline with custom data
To use your own data, place your `.tar` archive of BAM files (aligned to HG19) in the `assets` directory.
Update `input_data_path` and `outdir` in `nextflow.config` to match your input file and desired output location.


## Citation
The methodology and probabilistic model for ichorCNA are described in:
Adalsteinsson, Ha, Freeman, et al. Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. (2017) Nature Communications Nov 6;8(1):1324. doi: 10.1038/s41467-017-00965-y

[ichorCNA GitHub repository](https://github.com/broadinstitute/ichorCNA)