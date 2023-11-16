# Deep mutational scanning of SARS-CoV-2 XBB.1.5 spike 
See [Dadonaite et al, bioRxiv, DOI 10.1101/2023.11.13.566961 (2023)](https://doi.org/10.1101/2023.11.13.566961
) for the paper describing this study.

For documentation of the analysis, see [https://dms-vep.github.io/SARS-CoV-2_XBB.1.5_spike_DMS/](https://dms-vep.github.io/SARS-CoV-2_XBB.1.5_spike_DMS/)

## Organization of this repo

### `dms-vep-pipeline-3` submodule

Most of the analysis is done by the [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), which was added as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to this pipeline via:

    git submodule add https://github.com/dms-vep/dms-vep-pipeline-3

This added the file [.gitmodules](.gitmodules) and the submodule [dms-vep-pipeline-3](dms-vep-pipeline-3), which was then committed to the repo.
Note that if you want a specific commit or tag of [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) or to update to a new commit, follow the [steps here](https://stackoverflow.com/a/10916398), basically:

    cd dms-vep-pipeline-3
    git checkout <commit>

and then `cd ../` back to the top-level directory, and add and commit the updated `dms-vep-pipeline-3` submodule.
You can also make changes to the [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) that you commit back to that repo.

### Additional steps specific to this pipeline
In addition to the main [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), there are some custom `snakemake` rules in [custom_rules.smk](custom_rules.smk) and the associated notebooks in [./notebooks/](notebooks) that run some rules and analyses specific to this repo.

### Multinomial logistic clade growth data
An analysis by Ben Murrell of SARS-CoV-2 clade growth using multinomial logistic models in the repo at [https://github.com/MurrellGroup/MultinomialLogisticGrowth](https://github.com/MurrellGroup/MultinomialLogisticGrowth) is added as a submodule in this repo, and is used in some of the analyses.

### Configuration and running the pipeline
The configuration for the pipeline is in [config.yaml](config.yaml) and the files in [./data/](data) referenced therein.
To run the pipeline, do:

    snakemake -j 8 --use-conda -s dms-vep-pipeline-3/Snakefile

To run on the Hutch cluster via [slurm](https://slurm.schedmd.com/), you can run the file [run_Hutch_cluster.bash](run_Hutch_cluster.bash):

    sbatch -c 8 run_Hutch_cluster.bash

### Results and documentation
The results of running the pipeline are placed in [./results/](results).
Only some of these results are tracked to save space (see [.gitignore](.gitignore)).

The pipeline builds HTML documentation for the pipeline in [./docs/](docs), which is rendered via GitHub Pages at [https://dms-vep.github.io/SARS-CoV-2_XBB.1.5_spike_DMS/](https://dms-vep.github.io/SARS-CoV-2_XBB.1.5_spike_DMS/).

### Library design
The design of the mutant library is contained in [./library_design/](library_design).
That design is not part of the pipeline but contains code that must be run separately with its own [conda](https://docs.conda.io/) environment.
