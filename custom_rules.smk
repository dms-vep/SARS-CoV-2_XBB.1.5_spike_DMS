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


rule escape_at_key_sites:
    """Analyze and make logo plots of escape at key sites."""
    input:
        dms_csv="results/summaries/summary.csv",
        per_antibody_csv="results/summaries/per_antibody_escape.csv",
        nb="notebooks/escape_at_key_sites.ipynb",
    output:
        nb="results/notebooks/escape_at_key_sites.ipynb",
    params:
        yaml=lambda _, input: yaml.round_trip_dump(
            {
                "pango_consensus_seqs_json": "https://raw.githubusercontent.com/corneliusroemer/pango-sequences/c64ef05e53debaa9cc65dd56d6eb83e31517179c/data/pango-consensus-sequences_summary.json",
                "dms_csv": input.dms_csv,
                "per_antibody_csv": input.per_antibody_csv,
            }
        ),
    log:
        log="results/logs/escape_at_key_sites.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        "papermill {input.nb} {output.nb} -y '{params.yaml}' &>> {log}"


rule compare_binding:
    """Compare ACE2 binding across datasets."""
    input:
        xbb_spike_csv="results/summaries/summary.csv",
        nb="notebooks/compare_binding.ipynb",
    params:
        yaml=lambda wc, input: yaml.round_trip_dump(
            {
                # ----------------------------------------
                # parameters for plots
                # ----------------------------------------
                "init_min_func_effect": -1.5,
                "clip_binding_upper": 4,
                "clip_binding_lower": -6,
                # ----------------------------------------
                # Other deep mutational scanning datasets
                # ----------------------------------------
                # Tyler Starr yeast display RBD DMS
                "starr_rbd_affinity":
                    "https://media.githubusercontent.com/media/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ/main/results/final_variant_scores/final_variant_scores.csv",
                # BA.2 in full spike DMS
                "ba2_spike_csv":
                    "https://raw.githubusercontent.com/dms-vep/SARS-CoV-2_Omicron_BA.2_spike_ACE2_binding/main/results/summaries/summary.csv",
                # XBB.1.5 in RBD DMS in lentiviral system
                "xbb_rbd_csv":
                    "https://raw.githubusercontent.com/dms-vep/SARS-CoV-2_XBB.1.5_RBD_DMS/main/results/summaries/summary.csv",
            }
            | {key: val for (key, val) in dict(input).items() if key != "nb"}
        ),
    output:
        merged_binding_csv="results/binding_comparison/merged_binding.csv",
        nb="results/notebooks/compare_binding.ipynb",
        binding_corr="results/binding_comparison/binding_corr.html",
        binding_dist="results/binding_comparison/binding_dist.html",
        binding_entry_corr="results/binding_comparison/binding_entry_corr.html",
        binding_escape_corr="results/binding_comparison/binding_ecape_corr.html",
    log:
        log="results/logs/compare_binding.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        """
        papermill {input.nb} {output.nb} \
            -y '{params.yaml}' \
            -p merged_binding_csv {output.merged_binding_csv} \
            -p binding_corr_html {output.binding_corr} \
            -p binding_dist_html {output.binding_dist} \
            -p binding_entry_corr_html {output.binding_entry_corr} \
            -p binding_escape_corr_html {output.binding_escape_corr} \
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
                "max_effect_std": 1.6,
                "init_min_times_seen": 4,
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
                "init_min_func_effect": -1.5,
                "max_effect_std": 1.6,
                "init_min_times_seen": 4,
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
        growth_rates_csv="MultinomialLogisticGrowth/model_fits/rates.csv",
    output:
        nb="results/notebooks/compare_natural.ipynb",
        pango_consensus_seqs_json="results/compare_natural/pango-consensus-sequences_summary.json",
        pair_growth_dms_csv="results/compare_natural/clade_pair_growth_dms.csv",
        clade_growth_dms_csv="results/compare_natural/clade_growth_dms.csv",
    params:
        pango_consensus_seqs_json="https://raw.githubusercontent.com/corneliusroemer/pango-sequences/c64ef05e53debaa9cc65dd56d6eb83e31517179c/data/pango-consensus-sequences_summary.json",
        yaml=lambda _, input, output: yaml.round_trip_dump(
            {
                "starting_clades": ["XBB"],  # clades descended from this
                "muts_to_toggle": ["L455F"],  # epistasis in affinity for L455F in FLiP
                "min_sequences": 400,  # require this many sequences per clade to use
                "split_by_rbd": False,  # whether to treat RBD and non-RBD mutations separately
                "dms_clade": "XBB.1.5",  # clade used for DMS
                "pair_min_spike_muts": 1,  # require clade pairs to have >= this many spike mutations
                "pair_max_spike_muts": None,  # require clade pairs to have <= this many spike mutations
                "n_random": 100,  # compute P values with this many randomizations of DMS data
                "phenotype_basic_colors": {
                    # define phenotypes and their colors. "basic" means not split by RBD, which is done
                    # later in code if `split_by_rbd`
                    "sera escape": "red",
                    "ACE2 binding": "blue",
                    "cell entry": "purple",
                },
                "exclude_clades": [],
                "growth_rates_csv": input.growth_rates_csv,
                "dms_summary_csv": input.dms_summary_csv,
                "pango_consensus_seqs_json": output.pango_consensus_seqs_json,
                "pair_growth_dms_csv": output.pair_growth_dms_csv,
                "clade_growth_dms_csv": output.clade_growth_dms_csv,
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

rule non_rbd_binding_natural:
    """Look at non-RBD mutation effects on ACE2 binding in natural viruses."""
    input:
        dms_summary_csv="results/summaries/summary.csv",
        pango_consensus_seqs_json=rules.compare_natural.output.pango_consensus_seqs_json,
        nb="notebooks/non_rbd_binding_natural.ipynb",
    output:
        nb="results/notebooks/non_rbd_binding_natural.ipynb",
    params:
        yaml=lambda wc, input: yaml.round_trip_dump(
            {
                "pango_consensus_seqs_json": input.pango_consensus_seqs_json,
                "xbb15_dms_csv": input.dms_summary_csv,
                "ba2_dms_csv":
                    "https://raw.githubusercontent.com/dms-vep/SARS-CoV-2_Omicron_BA.2_spike_ACE2_binding/main/results/summaries/summary.csv",
            }
        ),
    log:
        log="results/logs/non_rbd_binding_natural.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        "papermill {input.nb} {output.nb} -y '{params.yaml}' &>> {log}"


rule func_effects_dist:
    """Distribution of functional effects and correlation with natural sequences."""
    input:
        xbb15_func_effects_csv="results/func_effects/averages/293T_high_ACE2_entry_func_effects.csv",
        site_numbering_map_csv=config["site_numbering_map"],
        nb="notebooks/func_effects_dist.ipynb",
    output:
        strain_corr="results/func_effects_analyses/strain_corr.html",
        natural_corr="results/func_effects_analyses/natural_corr.html",
        effects_boxplot="results/func_effects_analyses/effects_boxplot.html",
        key_muts_plot="results/func_effects_analyses/key_mutations.html",
        nb="results/notebooks/func_effects_dist.ipynb",
    params:
        yaml=lambda _, input, output: yaml.round_trip_dump(
            {
                "fitness_csv": "https://raw.githubusercontent.com/jbloomlab/SARS2-mut-fitness/main/results_public_2023-10-01/aa_fitness/aa_fitness.csv",
                "xbb15_func_effects_csv": input.xbb15_func_effects_csv,
                "ba2_func_effects_csv":
                    "https://raw.githubusercontent.com/dms-vep/SARS-CoV-2_Omicron_BA.2_spike_ACE2_binding/main/results/func_effects/averages/293T_high_ACE2_entry_func_effects.csv",
                "site_numbering_map_csv": input.site_numbering_map_csv,
                "init_min_times_seen": 3,
                "init_min_n_libraries": 2,
                "init_expected_count": 20,
                "max_effect_std": 1.6,
                "key_mutations": ["P1143L", "F456L", "V483-"],
                "strain_corr_html": output.strain_corr,
                "natural_corr_html": output.natural_corr,
                "effects_boxplot_html": output.effects_boxplot,
                "key_muts_html": output.key_muts_plot,
            }
        ),
    log:
        "results/logs/func_effects_dist.txt",
    shell:
        "papermill {input.nb} {output.nb} -y '{params.yaml}' &>> {log}"


# Files (Jupyter notebooks, HTML plots, or CSVs) that you want included in
# the HTML docs should be added to the nested dict `docs`:
docs["Additional files and charts"] = {
    "Analysis of escape and other properties at key sites": {
        "Notebook performing analysis and making logo plots": rules.escape_at_key_sites.output.nb,
    },
    "Analysis of ACE2 binding data and comparison to other experiments": {
        "Interactive charts": {
            "Correlations among experiments":
                rules.compare_binding.output.binding_corr,
            "Distribution of RBD and non-RBD ACE2 binding":
                rules.compare_binding.output.binding_dist,
            "Correlation of ACE2 binding to viral entry":
                rules.compare_binding.output.binding_entry_corr,
            "Correlation of ACE2 binding to viral escape":
                rules.compare_binding.output.binding_escape_corr,
        },
        "CSV of ACE2 binding from different experiments":
            rules.compare_binding.output.merged_binding_csv,
    },
    "Comparison of escape in medium and high ACE2 cells": {
        "Interactive chart comparing escape":
            rules.compare_high_medium_ace2_escape.output.chart,
    },
    "Comparison of escape in full-spike and RBD deep mutational scans": {
        "Interactive chart comparing escape in spike vs RBD scans":
            rules.compare_spike_rbd_escape.output.corr_chart,
        "Distributions of escape by RBD and non-RBD mutations in spike scan":
            rules.compare_spike_rbd_escape.output.dist_chart,
    },
    "DMS measurements versus clade growth": {
        "Notebook comparing change in clade growth to change in DMS phenotype":
            rules.compare_natural.output.nb,
        "CSV with data for comparison of changes in growth vs DMS phenotype for clade pairs":
            rules.compare_natural.output.pair_growth_dms_csv,
        "CSV with DMS phenotypes and growth for all individual clades":
            rules.compare_natural.output.clade_growth_dms_csv,
    },
    "Analysis of mutational effects on cell entry": {
        "Correlation of cell entry effects among strains": rules.func_effects_dist.output.strain_corr,
        "Correlation with fitness effects estimated from natural sequences":
            rules.func_effects_dist.output.natural_corr,
        "Distribution of cell entry effects": rules.func_effects_dist.output.effects_boxplot,
        "Effects of key mutations on cell entry": rules.func_effects_dist.output.key_muts_plot,
    },
    "ACE2 binding effects of non-RBD mutations in natural sequences": {
        "Notebook analyzing ACE2 binding effects of non-RBD mutations in natural sequences":
            rules.non_rbd_binding_natural.output.nb,
    },
    "Spike site numbering": {
        "CSV converting sequential sites in XBB.1.5 spike to Wuhan-Hu-1 reference sites":
            config["site_numbering_map"],
    },
}
