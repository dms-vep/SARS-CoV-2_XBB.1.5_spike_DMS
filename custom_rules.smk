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


rule binding_vs_escape:
    """Compare binding and escape at key sites."""
    input:
        dms_csv="results/summaries/summary.csv",
        nb="notebooks/binding_vs_escape.ipynb",
    output:
        nb="results/notebooks/binding_vs_escape.ipynb",
        logoplot_subdir=directory("results/binding_vs_escape/logoplots"),
    params:
        yaml=lambda _, input, output: yaml.round_trip_dump(
            {
                "dms_csv": input.dms_csv,
                "logoplot_subdir": output.logoplot_subdir,
            }
        ),
    log:
        log="results/logs/binding_vs_escape.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        "papermill {input.nb} {output.nb} -y '{params.yaml}' &> {log}"


rule escape_at_key_sites:
    """Analyze and make logo plots of escape at key sites."""
    input:
        dms_csv="results/summaries/summary.csv",
        per_antibody_csv="results/summaries/per_antibody_escape.csv",
        codon_seq="data/XBB_1_5_spike_codon.fa",
        nb="notebooks/escape_at_key_sites.ipynb",
    output:
        nb="results/notebooks/escape_at_key_sites.ipynb",
        logoplot_subdir=directory("results/key_sites/logoplots"),
    params:
        yaml=lambda _, input, output: yaml.round_trip_dump(
            {
                "dms_csv": input.dms_csv,
                "per_antibody_csv": input.per_antibody_csv,
                "codon_seq": input.codon_seq,
                "logoplot_subdir": output.logoplot_subdir,
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


rule pango_consensus_seqs_json:
    """Get JSON with pango consensus seqs."""
    params:
        pango_consensus_seqs_json="https://raw.githubusercontent.com/corneliusroemer/pango-sequences/9ef44f19bcea322c579b91e59756c4a27e7f943d/data/pango-consensus-sequences_summary.json",
    output:
        json="results/compare_natural/pango-consensus-sequences_summary.json",
    log:
        "results/logs/pango_consensus_seqs_json.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        "curl {params.pango_consensus_seqs_json} -o {output.json} &> {log}"


rule muts_per_clade:
    """Plot mutations per clade."""
    input:
        json=rules.pango_consensus_seqs_json.output.json,
        nb="notebooks/muts_per_clade.ipynb",
    output:
        nb="results/notebooks/muts_per_clade.ipynb",
        chart="results/compare_natural/muts_per_clade.html",
    log:
        "results/logs/muts_per_clade.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p chart_html {output.chart} \
            -p pango_consensus_seqs_json {input.json} \
            &> {log}
        """


# sets of measurements to compare to natural evolution
phenos_compare_natural = {
    "current_dms": {
        "input_data": "results/summaries/summary.csv",
        "rename_cols": {
            "human sera escape": "sera escape",
            "spike mediated entry": "cell entry",
        },
        "phenotype_colors": {
            "sera escape": "red",
            "ACE2 binding": "blue",
            "cell entry": "purple",
        },
        "title": "XBB.1.5 full-spike DMS phenotypes",
        "missing_muts": "drop",  # drop clades with missing muts
    },
    "yeast_RBD_DMS": {
        "input_data": "data/compare_natural_datasets/yeast_RBD_DMS.csv",
        "rename_cols": {},
        "phenotype_colors": {"escape": "red", "ACE2 affinity": "blue", "RBD expression": "purple"},
        "title": "yeast RBD DMS phenotypes",
        "missing_muts": "zero",  # set missing (non-RBD) mutations to zero
    },
    "muts_from_Wuhan-Hu-1": {
        "input_data": "data/compare_natural_datasets/incremental_Hamming_distance_from_Wuhan-Hu-1.csv",
        "rename_cols": {"incremental Hamming distance": "distance"},
        "phenotype_colors": {"distance": "gray"},
        "title": "relative distance from Wuhan-Hu-1",
        "missing_muts": "drop",  # drop clades with missing muts
    },
    "EVEscape": {
        "input_data": "data/compare_natural_datasets/EVEscape_XBB_single_mutation_predictions.csv",
        "rename_cols": {},
        "phenotype_colors": {"EVEscape": "gray"},
        "phenotype_colors": {"EVEscape": "gray"},
        "title": "EVEscape",
        "missing_muts": "drop",  # drop clades with missing muts
    },
    "EVEscape_components": {
        "input_data": "data/compare_natural_datasets/EVEscape_XBB_single_mutation_predictions.csv",
        "rename_cols": {"fitness_evol_indices": "EVE fitness", "dissimilarity_charge_hydrophobicity": "aa dissimilarity", "accessibility_wcn": "accessibility"},
        "phenotype_colors": {"EVE fitness": "red", "aa dissimilarity": "blue", "accessibility": "green"},
        "title": "EVEscape components",
        "missing_muts": "drop",  # drop clades with missing muts
    },
}

rule compare_natural:
    """Compare DMS (or other) phenotype measurements to natural sequence evolution."""
    input:
        input_data=lambda wc: phenos_compare_natural[wc.pheno]["input_data"],
        nb="notebooks/compare_natural.ipynb",
        growth_rates_csv="MultinomialLogisticGrowth/model_fits/rates.csv",
        pango_consensus_seqs_json=rules.pango_consensus_seqs_json.output.json,
    output:
        nb="results/notebooks/{pheno}_compare_natural.ipynb",
        pair_growth_dms_csv="results/compare_natural/{pheno}_clade_pair_growth.csv",
        clade_growth_dms_csv="results/compare_natural/{pheno}_clade_growth.csv",
        pair_corr_html="results/compare_natural/{pheno}_clade_pair_growth.html",
        clade_corr_html="results/compare_natural/{pheno}_clade_growth.html",
        pair_ols_html="results/compare_natural/{pheno}_ols_clade_pair_growth.html",
    params:
        yaml=lambda wc, input, output: yaml.round_trip_dump(
            {
                "starting_clades": ["XBB"],  # clades descended from this
                "exclude_muts": [],  # exclude clades w these mutations
                "min_sequences": 400,  # require this many sequences per clade to use
                "split_by_rbd": False,  # whether to treat RBD and non-RBD mutations separately
                "dms_clade": "XBB.1.5",  # clade used for DMS
                "pair_min_spike_muts": 1,  # require clade pairs to have >= this many spike mutations
                "pair_max_spike_muts": None,  # require clade pairs to have <= this many spike mutations
                "n_random": 100,  # compute P values with this many randomizations of DMS data
                # rename columns in input data
                "rename_cols": phenos_compare_natural[wc.pheno]["rename_cols"],
                "title": phenos_compare_natural[wc.pheno]["title"],
                # "basic" means not split by RBD, which is done later in code if `split_by_rbd`
                "phenotype_basic_colors": phenos_compare_natural[wc.pheno]["phenotype_colors"],
                # "drop" clades with missing mutations, or set missing mutations to "zero"
                "missing_muts": phenos_compare_natural[wc.pheno]["missing_muts"],
                "exclude_clades": [],  # exclude these clades
                "growth_rates_csv": input.growth_rates_csv,
                "input_data": input.input_data,
                "pango_consensus_seqs_json": input.pango_consensus_seqs_json,
                "pair_growth_dms_csv": output.pair_growth_dms_csv,
                "clade_growth_dms_csv": output.clade_growth_dms_csv,
                "pair_corr_html": output.pair_corr_html,
                "clade_corr_html": output.clade_corr_html,
                "pair_ols_html": output.pair_ols_html,
            }
        ),
    log:
        log="results/logs/{pheno}_compare_natural.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        "papermill {input.nb} {output.nb} -y '{params.yaml}' &> {log}"


rule compare_natural_ba2_ba5_xbb:
    """Compare DMS (or other) phenotype measurements to natural sequence evolution.
    
    Differs from `compare_natural` by using all BA.2, BA.5, and XBB-descended clades.
    But written as separate rule to preserve backward compatibility with file-naming
    used in repo for original version.

    """
    input:
        input_data=lambda wc: phenos_compare_natural[wc.pheno]["input_data"],
        nb="notebooks/compare_natural.ipynb",
        growth_rates_csv="MultinomialLogisticGrowth/model_fits/rates.csv",
        pango_consensus_seqs_json=rules.pango_consensus_seqs_json.output.json,
    output:
        nb="results/notebooks/{pheno}_compare_natural_ba2_ba5_xbb.ipynb",
        pair_growth_dms_csv="results/compare_natural/{pheno}_clade_pair_growth_ba2_ba5_xbb.csv",
        clade_growth_dms_csv="results/compare_natural/{pheno}_clade_growth_ba2_ba5_xbb.csv",
        pair_corr_html="results/compare_natural/{pheno}_clade_pair_growth_ba2_ba5_xbb.html",
        clade_corr_html="results/compare_natural/{pheno}_clade_growth_ba2_ba5_xbb.html",
        pair_ols_html="results/compare_natural/{pheno}_ols_clade_pair_growth_ba2_ba5_xbb.html",
    params:
        yaml=lambda wc, input, output: yaml.round_trip_dump(
            {
                "starting_clades": ["BA.2", "BA.5", "XBB"],  # clades descended from this
                "exclude_muts": [],  # exclude clades w these mutations
                "min_sequences": 400,  # require this many sequences per clade to use
                "split_by_rbd": False,  # whether to treat RBD and non-RBD mutations separately
                "dms_clade": "XBB.1.5",  # clade used for DMS
                "pair_min_spike_muts": 1,  # require clade pairs to have >= this many spike mutations
                "pair_max_spike_muts": None,  # require clade pairs to have <= this many spike mutations
                "n_random": 100,  # compute P values with this many randomizations of DMS data
                # rename columns in input data
                "rename_cols": phenos_compare_natural[wc.pheno]["rename_cols"],
                "title": phenos_compare_natural[wc.pheno]["title"],
                # "basic" means not split by RBD, which is done later in code if `split_by_rbd`
                "phenotype_basic_colors": phenos_compare_natural[wc.pheno]["phenotype_colors"],
                # "drop" clades with missing mutations, or set missing mutations to "zero"
                "missing_muts": phenos_compare_natural[wc.pheno]["missing_muts"],
                "exclude_clades": [],  # exclude these clades
                "growth_rates_csv": input.growth_rates_csv,
                "input_data": input.input_data,
                "pango_consensus_seqs_json": input.pango_consensus_seqs_json,
                "pair_growth_dms_csv": output.pair_growth_dms_csv,
                "clade_growth_dms_csv": output.clade_growth_dms_csv,
                "pair_corr_html": output.pair_corr_html,
                "clade_corr_html": output.clade_corr_html,
                "pair_ols_html": output.pair_ols_html,
            }
        ),
    log:
        log="results/logs/{pheno}_compare_natural_ba2_ba5_xbb.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        "papermill {input.nb} {output.nb} -y '{params.yaml}' &> {log}"


rule compare_ba_2_86:
    """Compare predicted phenotypes of actual and randomized sequences related to BA.2.86."""
    input:
        clade_phenotypes_csv="SARS2-spike-predictor-phenos/results/clade_phenotypes.csv",
        mutation_phenotypes_csv="SARS2-spike-predictor-phenos/results/mutation_phenotypes.csv",
        gisaid_mutation_counts_csv="data/GISAID_alignment_counts_2024-01-27.csv",
        nb="notebooks/compare_BA.2.86.ipynb",
    params:
        yaml=lambda _, input: yaml.round_trip_dump(
            {
                "gisaid_min_counts": 50,  # draw random mutations from those with >= this many GISAID counts
                "nrandom": 1000,  # number randomized sequences, for descendants it is 10x less
                "linear_models": {  # weights of phenotypes in linear model
                    "spike pseudovirus DMS (combined phenotypes)":
                        {
                            "spike pseudovirus DMS human sera escape": 38,
                            "spike pseudovirus DMS ACE2 binding": 2,
                            "spike pseudovirus DMS spike mediated entry": 16,
                        },
                    "RBD yeast-display DMS (combined phenotypes)":
                        {
                            "RBD yeast-display DMS ACE2 affinity": 19,
                            "RBD yeast-display DMS RBD expression": 20,
                            "RBD yeast-display DMS escape": 27,
                        },
                },
                "clade_phenotypes_csv": input.clade_phenotypes_csv,
                "mutation_phenotypes_csv": input.mutation_phenotypes_csv,
                "gisaid_mutation_counts_csv": input.gisaid_mutation_counts_csv,
            }
        ),
    output:
        nb="results/notebooks/compare_BA.2.86.ipynb",
    log:
        log="results/logs/compare_BA.2.86.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    shell:
        "papermill {input.nb} {output.nb} -y '{params.yaml}' &> {log}"


rule non_rbd_binding_natural:
    """Look at non-RBD mutation effects on ACE2 binding in natural viruses."""
    input:
        dms_summary_csv="results/summaries/summary.csv",
        pango_consensus_seqs_json=rules.pango_consensus_seqs_json.output.json,
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
        "papermill {input.nb} {output.nb} -y '{params.yaml}' &> {log}"


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


rule escape_by_prior_infections:
    """Escape stratified by prior infections."""
    # somewhat hackily copied from the pipeline summary plot
    input:
        **dict(rules.summary.input),
    output:
        chart_overlaid="results/escape_by_prior_infections/summary_overlaid_nolegend.html",
        chart_faceted="results/escape_by_prior_infections/summary_faceted_nolegend.html",
        csv="results/escape_by_prior_infections/summary.csv",
        per_antibody_escape_csv="results/escape_by_prior_infections/per_antibody_escape.csv",
        nb="results/notebooks/escape_by_prior_infections.ipynb",
    params:
        yaml=lambda _, input: yaml.round_trip_dump(
            {
                "config": (
                    {
                        key: val
                        for (key, val) in summary_config.items()
                        if key not in ["title", "legend", "antibody_escape"]
                    }
                    | {
                        "antibody_escape" : {
                            "one infection": {
                                "stat": "escape_median",
                                "positive_color": "#56B4E9",
                                "negative_color": "#E69F00",
                                "max_at_least": 1,
                                "min_at_least": -1,
                                "le_filters": {"escape_std": 1.5},
                                "antibody_list": {
                                    "sera_493C_mediumACE2": "serum 493C",
                                    "sera_498C_mediumACE2": "serum 498C",
                                    "sera_500C_mediumACE2": "serum 500C",
                                    "sera_501C_mediumACE2": "serum 501C",
                                    "sera_503C_mediumACE2": "serum 503C",
                                    "sera_505C_mediumACE2": "serum 505C",
                                },
                            },
                            "multiple infections": {
                                "stat": "escape_median",
                                "positive_color": "#56B4E9",
                                "negative_color": "#E69F00",
                                "max_at_least": 1,
                                "min_at_least": -1,
                                "le_filters": {"escape_std": 1.5},
                                "antibody_list": {
                                    "sera_287C_mediumACE2": "serum 287C",
                                    "sera_288C_mediumACE2": "serum 288C",
                                    "sera_343C_mediumACE2": "serum 343C",
                                    "sera_497C_mediumACE2": "serum 497C",
                                },
                            },
                        },
                    }
                ),
                "input_csvs": dict(input),
            }
        ),
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    log:
        "results/logs/escape_by_prior_infections.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p site_numbering_map_csv {input.site_numbering_map_csv} \
            -p chart_faceted {output.chart_faceted} \
            -p chart_overlaid {output.chart_overlaid} \
            -p output_csv_file {output.csv} \
            -p per_antibody_escape_csv {output.per_antibody_escape_csv} \
            -y "{params.yaml}" \
            &> {log}
        """

# Files (Jupyter notebooks, HTML plots, or CSVs) that you want included in
# the HTML docs should be added to the nested dict `docs`:
docs["Additional files and charts"] = {
    "Analysis of escape and other properties at key sites": {
        "Notebook making logo plots of escape at key sites": rules.escape_at_key_sites.output.nb,
        "Notebook comparing binding vs escape at key sites": rules.binding_vs_escape.output.nb,
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
    "DMS measurements versus clade growth (XBB clades)": {
        f"Comparison for {pheno} (XBB clades)": {
            "Correlation of change in clade growth versus phenotype for clade pairs":
                rules.compare_natural.output.pair_corr_html.format(pheno=pheno),
            "OLS change in clade growth versus phenotype for clade pairs":
                rules.compare_natural.output.pair_ols_html.format(pheno=pheno),
            "Correlation of absolute clade growth versus phenotype":
                rules.compare_natural.output.clade_corr_html.format(pheno=pheno),
            "Notebook comparing change in clade growth to change in phenotype":
                rules.compare_natural.output.nb.format(pheno=pheno),
            "CSV with data for comparison of changes in growth vs phenotype for clade pairs":
                rules.compare_natural.output.pair_growth_dms_csv.format(pheno=pheno),
            "CSV with phenotypes and growth for all individual clades":
                rules.compare_natural.output.clade_growth_dms_csv.format(pheno=pheno),
        }
        for pheno in phenos_compare_natural
    },
    "DMS measurements versus clade growth (BA.2, BA.5, and XBB clades)": {
        f"Comparison for {pheno} (BA.2, BA.5, and XBB clades)": {
            "Correlation of change in clade growth versus phenotype for clade pairs":
                rules.compare_natural_ba2_ba5_xbb.output.pair_corr_html.format(pheno=pheno),
            "OLS change in clade growth versus phenotype for clade pairs":
                rules.compare_natural_ba2_ba5_xbb.output.pair_ols_html.format(pheno=pheno),
            "Correlation of absolute clade growth versus phenotype":
                rules.compare_natural_ba2_ba5_xbb.output.clade_corr_html.format(pheno=pheno),
            "Notebook comparing change in clade growth to change in phenotype":
                rules.compare_natural_ba2_ba5_xbb.output.nb.format(pheno=pheno),
            "CSV with data for comparison of changes in growth vs phenotype for clade pairs":
                rules.compare_natural_ba2_ba5_xbb.output.pair_growth_dms_csv.format(pheno=pheno),
            "CSV with phenotypes and growth for all individual clades":
                rules.compare_natural_ba2_ba5_xbb.output.clade_growth_dms_csv.format(pheno=pheno),
        }
        for pheno in phenos_compare_natural
    },
    "Comparison to BA.2.86 evolution": {
        "Notebook comparing phenotypes to BA.2.86 evolution":
            rules.compare_ba_2_86.output.nb,
    },
    "Escape stratified by one or multiple prior infections": {
        "Notebook plotting escape with one versus multiple prior infections": 
            rules.escape_by_prior_infections.output.nb,
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
