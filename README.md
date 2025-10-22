# MagicWise Index

This tool is intended to be used to used to create a comprehensive score to compare two or more groups and identify causal genes. This is done in two parts with BeWISE (Burden Estimate from Weighted Integration of Site-specific Epigenetic Changes) and BeMAGIC (Burden Estimate from Modified and Associated Genetic Change) to create the MagicWise Index.

## Getting Started

### Installation

Combined Burden Estimate is run using [nextflow](https://www.nextflow.io). There are several ways to install and manage nextflow as noted on their website. Once installed, run the following to insure proper setup. 

```
# Run from directory that nextflow is installed in
./nextflow run hello
```

To download the code to run Combined Burden Estimate, run

```
git clone https://github.com/chloemoel/MagicWiseIndex

cd MagicWiseIndex
```
To install necessary data files, run: 

```
./install.sh
```

This workflow is written to work with a container system. There is no container install needed before running the pipeline, as Nextflow will find, install, and cache the container. The example config file mentioned below shows how to use and set a cache directory. This workflow will work with Apptainer, Singularity, or Docker systems.

### Config Files
This pipeline is run dependent on a nextflow config file saved in the base directory of the workflow. You can find an example nexflow config file [here](nextflow.config.example)

You can find more information on config files [here](https://www.nextflow.io/docs/latest/config.html)

### Input File formats

* sample_sheet (csv file) -- this sample sheet is used for the BeWISE calculation only, as the genetic score depends on a vcf with headers as study id names. This sample sheet is based on the Illumina sample sheet format with headers as follows: Sample_Name,Sample_Well,Sample_Plate,Sample_Group,Pool_ID,Sentrix_ID,Sentrix_Position. Only Sample_Name, Sentrix_ID, and Sentrix_Position are required. Sample name MUST match the headers of the vcf file. 
* additional_data (csv file) -- this a sheet with additional data used for batch correction during the BeWISE calculation. This csv file must contain the study id in the first column, and then any other information in subsequent columns
* sample_m_vals (csv file) -- a file with probes as rows and samples as columns. Samples should be in the [Sentrix_ID]_[Sentrix_Position] format. To calculate m values, we use `SeSAME` (more info can be found [here](https://zhou-lab.github.io/sesame/v1.16/sesame.html)). Probes should be in either the EPIC or 450k array format (NOT EPICv2). 
* vcf -- a multisample vcf file with all samples included in a single file. QC for genotyping calls should be done before running this tool.
* batch_correction -- list of variables to correct for in batch correction during BeWISE. Should be either in the sample_sheet file (ie, Sentrix_ID for chip correction) or a column from the additional data file. 

### Useage

To run this pipeline after you have installed nextflow, run the install script, and have created your config file:

```
./path/to/nextflow run main.nf -profile high_computing_cluster
```

If you are running after editing a file, you can run the following to resume

```
./path/to/nextflow run main.nf -profile high_computing_cluster -resume
```

Your score outputs can be found in the `outputs` directory

To run the example data, run either command (depending on your container system)
```
./path/to/nextflow -C nextflow.config.example run main.nf -profile docker 
./path/to/nextflow -C nextflow.config.example run main.nf -profile singularity
```
