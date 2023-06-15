# Deep mutational scanning of SARS-CoV-2 XBB.1.5 spike 
Study by Bernadeta Dadonaite and Jesse Bloom.

For documentation of the analysis, see [https://dms-vep.github.io/https://github.com/dms-vep/SARS-CoV-2_XBB.1.5_spike_DMS/](https://dms-vep.github.io/https://github.com/dms-vep/SARS-CoV-2_XBB.1.5_spike_DMS/).

## Organization of this repo

### `dms-vep-pipeline` submodule

Most of the analysis is done by the [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline), which was added as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to this pipeline via:

    git submodule add https://github.com/dms-vep/dms-vep-pipeline

This added the file [.gitmodules](.gitmodules) and the submodule [dms-vep-pipeline](dms-vep-pipeline), which was then committed to the repo.
Note that if you want a specific commit or tag of [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) or to update to a new commit, follow the [steps here](https://stackoverflow.com/a/10916398), basically:

    cd dms-vep-pipeline
    git checkout <commit>

and then `cd ../` back to the top-level directory, and add and commit the updated `dms-vep-pipeline` submodule.
You can also make changes to the [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) that you commit back to that repo.

### Code and configuration
The [snakemake](https://snakemake.readthedocs.io/) pipeline itself is run by the [Snakefile](Snakefile), which includes [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) reads its configuration from [config.yaml](config.yaml).
The [conda](https://docs.conda.io/) environment used by the pipeline is that specified in the `environment.yml` file in [dms-vep-pipeline](dms-vep-pipeline).

Additional scripts and notebooks that are specific to this analysis and not part of [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) are in [./scripts/](scripts) and [./notebooks/](notebooks).

### Input data
Input data for the pipeline are in [./data/](data).

### Results and documentation
The results of running the pipeline are placed in [./results/](results).
Only some of these results are tracked to save space (see [.gitignore](.gitignore)).

The pipeline builds HTML documentation for the pipeline in [./docs/](docs), which is rendered via GitHub Pages at [https://dms-vep.github.io/SARS-CoV-2_Omicron_BA.2_spike_DMS/](https://dms-vep.github.io/SARS-CoV-2_Omicron_BA.2_spike_DMS/).

### Library design
The design of the mutant library is contained in [./library_design/](library_design).
That design is not part of the pipeline but contains code that must be run separately with its own [conda](https://docs.conda.io/) environment.

## Running the pipeline
To run the pipeline, build the conda environment `dms-vep-pipeline` in the `environment.yml` file of [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline), activate it, and run [snakemake](https://snakemake.readthedocs.io/), such as:

    conda activate dms-vep-pipeline
    snakemake -j 32 --use-conda --rerun-incomplete

To run on the Hutch cluster via [slurm](https://slurm.schedmd.com/), you can run the file [run_Hutch_cluster.bash](run_Hutch_cluster.bash):

    sbatch -c 32 run_Hutch_cluster.bash
