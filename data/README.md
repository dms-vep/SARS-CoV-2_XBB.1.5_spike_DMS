# Input data

## PacBio full-length variant sequencing
[PacBio_runs.csv](PacBio_runs.csv) contains information on PacBio runs linking the barcodes to the mutations. It has the following columns:
 - `library`: name of the library sequenced
 - `run`: date of the pacbio library submission (use this date to refer to [experimental notebook](https://docs.google.com/document/d/1kWsdfg_-4m59Jj2jqEZwGwP5sN3cCqypFzaeGgUX144/edit?usp=sharing))
 - `fastq`: FASTQ file from running CCS

The amplicons itself is defined in [PacBio_amplicon.gb](PacBio_amplicon.gb), and [PacBio_feature_parse_specs.yaml](PacBio_feature_parse_specs.yaml) specifies how to parse that with [alignparse](https://jbloomlab.github.io/alignparse/).

## Illumina barcode sequencing
[barcode_runs.csv](barcode_runs.csv) has the Illumina barcode-sequencing runs used to count barcodes in different conditions.
It describes the samples, which should be named as clearly as possible.
It has the following columns:

 - *date*: date experiment was performed in `YYYY-MM-DD` encoding.
 - *virus_batch*: batch of virus used for the experiment.
 - *library*: which virus library was used.
 - *sample_type*: can be one of the following:
   + *VSVG_control*: entry mediated by VSVG
   + *no-antibody_control*: entry mediated by VEP of interest
   + *antibody*: encompasses sera and antibodies
 - *antibody*: name of the antibody if this sample has *sample_type* of *antibody*
 - *antibody_concentration*: concentration of antibody if this sample has *sample_type* of antibody. For sera, should be a fraction < 1 giving dilution (**not** a dilution factor).
 - *replicate*: experimental replicate.
 - *fastq_R1*: path to R1 FASTQ file, or semi-colon de-limited list of multiple FASTQs
 - *exclude_after_counts*: set to *yes* if barcode run should be excluded after counting barcodes
 - *notes*: any other notes about the sample.

The file [neutralization_standard_barcodes.csv](neutralization_standard_barcodes.csv) has the barcodes for the neutralization standard that is not expected to be neutralized by the sera or antibodies.

## Configuration for `polyclonal` analysis
[polyclonal_config.yaml](polyclonal_config.yaml) specifies how the analysis with [polyclonal](https://jbloomlab.github.io/polyclonal/) is done.
For each antibody listed in [barcode_runs.csv](barcode_runs.csv), specify:

 - *max_epitopes*: the maximum number of epitopes to test. The fitting keeps testing more epitopes up to this max until additional epitopes have no additional value.
 - *n_bootstrap_samples*: number of bootstrap samples to use.
 - *reg_escape_weight*: regularization weight for mutation-escape values.
 - *reg_spread_weight*: regularization weight for spread of escape values at each site.
 - *reg_activity_weight*: regularization weight for epitope activities.
 - *times_seen*: the `times_seen` value used for plotting the results (number of variants a mutation must be found in).
 - *min_epitope_activity_to_include*: keep adding epitopes until activity <= this.
