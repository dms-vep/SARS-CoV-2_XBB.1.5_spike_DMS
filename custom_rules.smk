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
                "init_min_func_effect": -1.5,
                "clip_affinity_upper": 4,
                "clip_affinity_lower": -6,
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


rule compare_high_medium_ace2_escape:
    """Compare escape on high- and medium-ACE2 cells."""
    input:
        escape=expand(
            rules.avg_escape.output.effect_csv,
            assay=["antibody_escape"],
            antibody=avg_assay_config["antibody_escape"],
        ),
        site_numbering_map=config["site_numbering_map"],
        func_effects="results/func_effects/averages/293T_high_ACE2_entry_func_effects.csv",
        nb="notebooks/compare_high_medium_ace2_escape.ipynb",
    output:
        chart="results/escape_comparisons/compare_high_medium_ace2_escape.html",
        nb="results/notebooks/compare_high_medium_ace2_escape.ipynb",
    params:
        yaml=lambda wc, input: yaml.round_trip_dump(
            {
                "init_min_func_effect": -2,
                "init_min_times_seen": 3,
                "init_floor_at_zero": False,
                "init_site_escape_stat": "sum",
                "escape_csvs": list(input.escape),
            }
        ),
    log:
        log="results/logs/compare_high_medium_ace2_escape.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        """
        papermill {input.nb} {output.nb} \
            -y '{params.yaml}' \
            -p site_numbering_map_csv {input.site_numbering_map} \
            -p func_effects_csv {input.func_effects} \
            -p chart_html {output.chart} \
            &> {log}
        """


rule compare_spike_rbd_escape:
    """Compare escape from full spike and RBD libraries."""
    input:
        spike_escape=expand(
            rules.avg_escape.output.effect_csv,
            assay=["antibody_escape"],
            antibody=avg_assay_config["antibody_escape"],
        ),
        site_numbering_map=config["site_numbering_map"],
        func_effects="results/func_effects/averages/293T_high_ACE2_entry_func_effects.csv",
        nb="notebooks/compare_spike_rbd_escape.ipynb",
    output:
        corr_chart="results/escape_comparisons/compare_spike_rbd_escape.html",
        dist_chart="results/escape_comparisons/rbd_vs_non_rbd_escape.html",
        nb="results/notebooks/compare_spike_rbd_escape.ipynb",
    params:
        yaml=lambda wc, input: yaml.round_trip_dump(
            {
                "init_min_func_effect": -2,
                "init_min_times_seen": 3,
                "init_floor_at_zero": False,
                "init_site_escape_stat": "sum",
                "spike_escape_csvs": list(input.spike_escape),
                "rbd_escape_csvs": [
                    os.path.join(
                        "https://raw.githubusercontent.com/dms-vep/SARS-CoV-2_XBB.1.5_RBD_DMS",
                        "main/results/antibody_escape/averages",
                        f"sera_{serum}_mediumACE2_mut_effect.csv",
                    )
                    for serum in ["493C", "498C", "500C", "503C", "343C"]
                ],
            }
        ),
    log:
        log="results/logs/compare_spike_rbd_escape.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        """
        papermill {input.nb} {output.nb} \
            -y '{params.yaml}' \
            -p site_numbering_map_csv {input.site_numbering_map} \
            -p func_effects_csv {input.func_effects} \
            -p corr_chart_html {output.corr_chart} \
            -p dist_chart_html {output.dist_chart} \
            &> {log}
        """


rule compare_natural:
    """Compare DMS measurements to natural sequence evolution."""
    input:
        dms_summary_csv="results/summaries/summary.csv",
        nb="notebooks/compare_natural.ipynb",
        growth_rates_csv="data/2023-09-18_Murrell_growth_estimates.csv",
    output:
        nb="results/notebooks/compare_natural.ipynb",
        pango_consensus_seqs_json="results/compare_natural/pango-consensus-sequences_summary.json",
        pango_dms_phenotypes_csv="results/compare_natural/pango_dms_phenotypes.csv",
        pango_by_date_html="results/compare_natural/pango_dms_phenotypes_by_date.html",
        pango_affinity_vs_escape_html="results/compare_natural/pango_affinity_vs_escape.html",
        pango_dms_vs_growth_regression_html="results/compare_natural/pango_dms_vs_growth_regression.html",
        pango_dms_vs_growth_regression_by_domain_html="results/compare_natural/pango_dms_vs_growth_regression_by_domain.html",
        pango_dms_vs_growth_corr_html="results/compare_natural/pango_dms_vs_growth_corr.html",
        pango_dms_vs_growth_corr_by_domain_html="results/compare_natural/pango_dms_vs_growth_corr_by_domain.html",
    params:
        pango_consensus_seqs_json="https://raw.githubusercontent.com/corneliusroemer/pango-sequences/67b13630832e163b9c0486dcce861d0106eaf07a/data/pango-consensus-sequences_summary.json",
        yaml=lambda _, input, output: yaml.round_trip_dump(
            {
                "starting_clades": ["BA.2", "BA.5", "XBB"],
                "dms_clade": "XBB.1.5",
                "dms_summary_csv": input.dms_summary_csv,
                "growth_rates_csv": input.growth_rates_csv,
                "pango_consensus_seqs_json": output.pango_consensus_seqs_json,
                "pango_dms_phenotypes_csv": output.pango_dms_phenotypes_csv,
                "pango_by_date_html": output.pango_by_date_html,
                "pango_affinity_vs_escape_html": output.pango_affinity_vs_escape_html,
                "pango_dms_vs_growth_regression_html": output.pango_dms_vs_growth_regression_html,
                "pango_dms_vs_growth_regression_by_domain_html": output.pango_dms_vs_growth_regression_by_domain_html,
                "pango_dms_vs_growth_corr_html": output.pango_dms_vs_growth_corr_html,
                "pango_dms_vs_growth_corr_by_domain_html": output.pango_dms_vs_growth_corr_by_domain_html,
                "n_random": 200,
                "exclude_clades": [],
                "exclude_clades_with_muts": [],
            }
        ),
    log:
        log="results/logs/compare_natural.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        """
        curl {params.pango_consensus_seqs_json} -o {output.pango_consensus_seqs_json} &> {log}
        papermill {input.nb} {output.nb} -y '{params.yaml}' &>> {log}
        """

rule non_rbd_affinity_natural:
    """Look at non-RBD mutation affects on ACE2 affinity in natural viruses."""
    input:
        dms_summary_csv="results/summaries/summary.csv",
        pango_consensus_seqs_json=rules.compare_natural.output.pango_consensus_seqs_json,
        nb="notebooks/non_rbd_affinity_natural.ipynb",
    output:
        nb="results/notebooks/non_rbd_affinity_natural.ipynb",
    params:
        yaml=lambda wc, input: yaml.round_trip_dump(
            {
                "pango_consensus_seqs_json": input.pango_consensus_seqs_json,
                "xbb15_dms_csv": input.dms_summary_csv,
                "ba2_dms_csv":
                    "https://raw.githubusercontent.com/dms-vep/SARS-CoV-2_Omicron_BA.2_spike_ACE2_affinity/main/results/summaries/summary.csv",
            }
        ),
    log:
        log="results/logs/non_rbd_affinity_natural.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        "papermill {input.nb} {output.nb} -y '{params.yaml}' &>> {log}"


# Files (Jupyter notebooks, HTML plots, or CSVs) that you want included in
# the HTML docs should be added to the nested dict `docs`:
docs["Additional files and charts"] = {
    "Analysis of ACE2 affinity data and comparison to other experiments": {
        "Interactive charts": {
            "Correlations among experiments":
                rules.compare_affinity.output.affinity_corr,
            "Distribution of RBD and non-RBD affinities":
                rules.compare_affinity.output.affinity_dist,
            "Correlation of affinity to viral entry":
                rules.compare_affinity.output.affinity_entry_corr,
            "Correlation of affinity to viral escape":
                rules.compare_affinity.output.affinity_escape_corr,
        },
        "CSV of affinities from different experiments":
            rules.compare_affinity.output.merged_affinity_csv,
    },
    "Comparison of escape in medium and high ACE2 cells": {
        "Interactive chart comparing escape":
            rules.compare_high_medium_ace2_escape.output.chart,
    },
    "Comparison of escape in full-spike and RBD deep mutational scans": {
        "Interactive chart comparing escape in spike vs RBD scans":
            rules.compare_spike_rbd_escape.output.corr_chart,
        "Distributions of escape by RBD and non-RBD mutations in spike scane":
            rules.compare_spike_rbd_escape.output.dist_chart,
    },
    "DMS phenotypes of natural Pango clades": {
        "Interactive chart of simple correlation of full-spike DMS phenotypes versus clade growth":
            rules.compare_natural.output.pango_dms_vs_growth_corr_html,
        "Interactive chart of linear regression of full-spike DMS phenotypes versus clade growth":
            rules.compare_natural.output.pango_dms_vs_growth_regression_html,
        "Interactive chart of simple correlation of RBD and non-RBD DMS phenotypes versus clade growth":
            rules.compare_natural.output.pango_dms_vs_growth_corr_by_domain_html,
        "Interactive chart of linear regression of RBD and non-RBD DMS phenotypes versus clade growth":
            rules.compare_natural.output.pango_dms_vs_growth_regression_by_domain_html,
        "Notebook analyzing DMS of natural Pango clades":
            rules.compare_natural.output.nb,
        "CSV with DMS phenotypes of the Pango clades":
            rules.compare_natural.output.pango_dms_phenotypes_csv,
    },
    "ACE2 affinity effects of non-RBD mutations in natural sequences": {
        "Notebook analyzing affinity effects of non-RBD mutations in natural sequences":
            rules.non_rbd_affinity_natural.output.nb,
    },
    "Spike site numbering": {
        "CSV converting sequential sites in XBB.1.5 spike to Wuhan-Hu-1 reference sites":
            config["site_numbering_map"],
    },
}
