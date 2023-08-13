# Input data
This subdirectory contains input data used by the pipeline.

## PacBio full-length variant sequencing to link barcodes

[PacBio_amplicon.gb](PacBio_amplicon.gb): Genbank file having features to parse with [alignparse](https://jbloomlab.github.io/alignparse/). Has H5 HA marked as *gene* (the gene of interest) and *barcode* site annotated.

[PacBio_feature_parse_specs.yaml](PacBio_feature_parse_specs.yaml): How to parse the PacBio amplicon using [alignparse](https://jbloomlab.github.io/alignparse/).

[PacBio_runs.csv](PacBio_runs.csv): List of PacBio CCS FASTQs used to link barcodes to variants.
It must have the following columns:

 - `library`: name of the library
 - `run`: name of the sequencing run, must be unique
 - `fastq`: FASTQ file from running CCS

## Site numbering
[site_numbering_map.csv](site_numbering_map.csv): Maps sequential 1, 2, ... numbering of the gene to a "reference" (WU-1) numbering scheme that represents the standard naming of sites for spike.
Also assigns each site to a regions of spike.
Has columns *sequential_site*, *reference_site*, and *region*, in this case *sequential_site* and *reference_site* are identical.

## Mutation-type classification
[data/mutation_design_classification.csv](data/mutation_design_classification.csv) classifies mutations into the different categories of designed mutations.
Has columns *sequential_site*, *amino_acid*, and *mutation_type*.

## Neutralization standard barcodes
[neutralization_standard_barcodes.csv](neutralization_standard_barcodes.csv) barcodes for the neutralization standards.
Has columns *barcode* and *name*, giving the barcode and name of this neutralization standard set.

## Barcode runs
[barcode_runs.csv](barcode_runs.csv) contains all samples and paths to sequencing files. It has the following format:

 - `sample`: sample name
 - `library`: name of library
 - `date`: date of sequencing
 - `fastq_R1`: path to one more FASTQ R1 sequencing files, multiple files should be semicolon-delimited

## Configuration for analyzing functional effects of mutations
[func_effects_config.yml](func_effects_config.yml) has the configuration for analyzing functional effects of mutations.
The format is explained within the file.

## Configuration for analyzing antibody escape
[antibody_escape_config.yml](antibody_escape_config.yml) has the configuration for analyzing effects of mutations on escape from antibodies or sera.
This same configuration can also be used to analyze escape from soluble receptor that inhibits entry to estimate how mutations affect receptor affinity.
The configuration for such receptor-inhibition experiments are also specified in this configuration as conceptually the experiments and analysis are basically the same, with just how the results are plotted and described different.
The format is explained within the file.
