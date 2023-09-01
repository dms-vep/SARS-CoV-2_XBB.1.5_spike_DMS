"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""


rule spatial_distances:
    """Get spatial distances from PDB."""
    input: 
        pdb="data/PDBs/aligned_spike_TM.pdb",
    output:
        csv="results/spatial_distances/spatial_distances.csv",
    params:
        target_chains=["A", "B", "C"],
    log:
        log="results/logs/spatial_distances.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    script:
        "scripts/spatial_distances.py"


rule compare_affinity:
    """Compare ACE2 affinity across datasets."""
    input:
        xbb_spike_csv="results/summaries/summary.csv",
        nb="notebooks/compare_affinity.ipynb",
    params:
        yaml=lambda wc, input: yaml.round_trip_dump(
            {
                # ----------------------------------------
                # parameters for plots
                # ----------------------------------------
                "init_min_func_effect": -2,
                "clip_affinity_upper": 2,
                "clip_affinity_lower": -4,
                # ----------------------------------------
                # Other deep mutational scanning datasets
                # ----------------------------------------
                # Tyler Starr yeast display RBD DMS
                "starr_rbd_affinity":
                    "https://media.githubusercontent.com/media/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ/main/results/final_variant_scores/final_variant_scores.csv",
                # BA.2 in full spike DMS
                "ba2_spike_csv":
                    "https://raw.githubusercontent.com/dms-vep/SARS-CoV-2_Omicron_BA.2_spike_ACE2_affinity/main/results/summaries/summary.csv",
                # XBB.1.5 in RBD DMS in lentiviral system
                "xbb_rbd_csv":
                    "https://raw.githubusercontent.com/dms-vep/SARS-CoV-2_XBB.1.5_RBD_DMS/main/results/summaries/summary.csv",
            }
            | {key: val for (key, val) in dict(input).items() if key != "nb"}
        ),
    output:
        merged_affinity_csv="results/affinity_comparison/merged_affinities.csv",
        nb="results/notebooks/compare_affinity.ipynb",
        affinity_corr="results/affinity_comparison/affinity_corr.html",
        affinity_dist="results/affinity_comparison/affinity_dist.html",
        affinity_entry_corr="results/affinity_comparison/affinity_entry_corr.html",
        affinity_escape_corr="results/affinity_comparison/affinity_ecape_corr.html",
    log:
        log="results/logs/compare_affinity.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        """
        papermill {input.nb} {output.nb} \
            -y '{params.yaml}' \
            -p merged_affinity_csv {output.merged_affinity_csv} \
            -p affinity_corr_html {output.affinity_corr} \
            -p affinity_dist_html {output.affinity_dist} \
            -p affinity_entry_corr_html {output.affinity_entry_corr} \
            -p affinity_escape_corr_html {output.affinity_escape_corr} \
            &> {log}
        """


# Files (Jupyter notebooks, HTML plots, or CSVs) that you want included in
# the HTML docs should be added to the nested dict `docs`:
docs["Additional files"] = {
    "Comparison of ACE2 affinities across experiments": {
        "Notebook comparing affinities": rules.compare_affinity.output.nb,
    },
    "Spike site numbering": {
        "CSV converting sequential sites in XBB.1.5 spike to Wuhan-Hu-1 reference sites":
            config["site_numbering_map"],
    },
}
