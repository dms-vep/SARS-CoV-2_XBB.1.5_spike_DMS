# Spike mutant library primer design

This repository designs primers for creating SARS-CoV-2 spike mutant libraries for deep mutational scanning.
The primers can be used for codon mutagenesis similar to that described [here](https://github.com/jbloomlab/CodonTilingPrimers), but is designed when wanting to create libraries that preferentially target mutationally tolerant and key sites similar in spirit to the approach described in [Dadonaite et al (2023)](https://www.sciencedirect.com/science/article/pii/S0092867423001034).

This repository does design for the XBB.1.5 spike.

Essentially, the pipeline identifies both specific mutations to target in the libraries, and sites that should undergo saturating mutagenesis via `NN[G/C]` primers.
Primers to make these mutations are then designed, potentially with different offsets.

## Design workflow

To run the workflow, build the `conda` environment in [environment.yml](environment.yml), then activate the file and run the `snakemake` pipeline in [Snakefile](Snakefile):

    conda activate spike_library_design
    snakemake -j 1
    
The results of running the workflow are then placed in [./results/](results).

The configuration for the workflow is in [config.yaml](config.yaml), input data are in [./data/](data), the scripts are in [./scripts/](scripts), and Jupyter notebooks are in [./notebooks/](notebooks).

The idea behind the design is as follows:

#### Specify the gene to mutagenize
We need the sequence of the gene to mutagenize.
Because we are designing primers that may extend outside the gene, this needs to be an "extended" sequence that also includes some upstream and downstream flanking.
A FASTA file containing the gene to mutagenize is specified via the `extended_gene` key in [config.yaml](config.yaml).
In that file, the gene sequence to be mutagenized should be in uppercase letters and flanking sequence not to be mutagenized should be in lower-case letters.

For instance, the gene that would be used for a XBB.1.5 full-spike library is in [data/extended_spike_XBB.1.5.fa](data/extended_spike_XBB.1.5.fa).

#### Get a mapping of sequential to reference site numbers
We want to refer to residues in reference (Wuhan-Hu-1) numbering.
This is done by aligning the protein we are mutagenizing with the Wuhan-Hu-1 reference specified under the `reference_gene` key in [config.yaml](config.yaml).
In this case, that file is [data/reference_spike.fa](data/reference_spike.fa).

The pipeline builds the site-numbering mapping, which is in [results/sequential_to_reference.csv](results/sequential_to_reference.csv).

#### Get GISAID alignment counts and UShER mutation counts
The pipeline gets the total counts of each non-reference amino-acid at each site in GISAID sequences from [here](https://mendel.bii.a-star.edu.sg/METHODS/corona/current/MUTATIONS/hCoV-19_Human_2019_WuhanWIV04/hcov19_Spike_mutations_table.html).

The pipeline also counts the number of unique occurrences on the tree (number of substitutions, not alignment counts) of each amino-acid mutation from the pre-built UShER tree.
These UShER mutation counts are done both for all SARS-CoV-2 clades, and for just the "recent" clades specified under `usher_recent_clades` in [config.yaml](config.yaml).

Information about the GISAId alignment counts, the overall UShER mutation counts, and the recent-clade UShER mutation counts is aggregated by the pipeline in the file [results/mutation_stats.csv](results/mutation_stats.csv).

#### Identify mutations to target and sites to saturate
The pipeline then identifies the mutations to target, as follows:

 1. Any mutation with counts >= any of the thresholds in the GISAID or UShER data listed under the `mutation_retain_thresholds` key in [config.yaml](config.yaml). However, deletion mutations are only included if they meet the counts **and** are in the range(s) specified under the `sites_to_allow_deletions` key in [config.yaml](config.yaml).
 2. Any mutation listed explicitly in `mutations_to_include` key in [config.yaml](config.yaml).
 3. All amino-acid and stop mutations at all of the sites to saturate. (Therefore, any mutation specified at a site to saturate is also a targeted mutation).
 
The pipeline writes the mutations to target to [results/targeted_mutations.csv](results/targeted_mutations.csv).
 
The sites to saturate are identified as follows:

  1. Any site with total mutation counts >= any of the thresholds specified under the `sites_to_allow_deletions` key in [config.yaml](config.yaml).
  2. All sites that differ between the sequence being mutagenized and the "reference" (eg, Wuhan-Hu-1 sequence) if the `saturate_diffs_from_reference` key in [config.yaml](config.yaml).
  3. All sites explicitly listed under the `sites_to_saturate` key in [config.yaml](config.yaml).
  
The pipeline writes the sites to saturate to [results/saturated_sites.csv](results/saturated_sites.csv).

The number of targeted mutations and saturated sites as a function of the parameters is shown in the interactive plot at [results/mutations_to_make.html](results/mutations_to_make.html).

#### Design the primers for oPools
The pipeline then designs actual primers that can be ordered in oPools to do the mutagenesis.

A total of four primer pools are designed:
 1. Forward primer pool for targeted mutations.
 2. Reverse primer pool for targeted mutations.
 3. Forward primer pool for saturated sites.
 4. Reverse primer pool for saturated sites.
 
For both the targeted and saturated sites, primers can be designed with one or more "offsets".
If the offset is zero, the mutation is placed at the center of the primer.
If the offset is negative, the mutation is shifted that many units towards the front of the primer; if it is positive the mutation is shifted towards the end of the primer.
The potential advantages of using offsets are (a) each mutation can then have slightly different primers which decreases chance of primer problems/failures, (b) offsets will alter how primers "overwrite" adjacent mutations and so increase chance of different combinations of nearby primers.

The offsets are specified under the `targeted_mutations_offsets` and `saturated_sites_offsets` keys in [config.yaml](config.yaml).
Primers are designed to have at least a melting temperature of `primer_min_tm`, a length of at least `primer_min_length`, and a length no more than `primer_max_length` as specified in [config.yaml](config.yaml).

The designed primer pools have a prefix in their name the string specified under `opool_prefix` in [config.yaml](config.yaml).

The designed pools themselves are in the file [results/oPools.csv](results/oPools.csv).

The number of unique sequences in each pool are in [results/oPool_stats.csv](results/oPool_stats.csv).
A recommended way to pool the primers is just to use concentrations of each pool proportional to the number of unique primers, so that each is at equal molarity.
This will give somewhat higher mutation rates at saturated sites since they are mutated by primers in both the targeted and saturated pools.
