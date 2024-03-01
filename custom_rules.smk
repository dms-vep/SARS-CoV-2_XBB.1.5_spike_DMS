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

# read the default SARS2-spike-predictor phenos YAML file and make replacements
sars2_spike_predictor_phenos_config_yaml = "SARS2-spike-predictor-phenos/config.yaml"
with open(sars2_spike_predictor_phenos_config_yaml) as f:
    sars2_spike_predictor_phenos_config = yaml.safe_load(f)
# specify specific versions of Pango clade data for reproducibility
sars2_spike_predictor_phenos_config["pango_json"] = "https://raw.githubusercontent.com/corneliusroemer/pango-sequences/71f1c4f249e2a6bc1d2a6b91a310ac1b610cd776/data/pango-consensus-sequences_summary.json"
sars2_spike_predictor_phenos_config["pango_growth_json"] = "https://data.nextstrain.org/files/workflows/forecasts-ncov/gisaid/pango_lineages/global/mlr/2024-02-27_results.json"
sars2_spike_predictor_phenos_outfiles = [
    "mutation_phenotypes_csv",
    "mutation_phenotypes_randomized_csv",
    "clade_phenotypes_csv",
    "clade_phenotypes_randomized_csv",
    "clade_phenotype_chart",
]
for sars2_spike_predictor_phenos_outfile in sars2_spike_predictor_phenos_outfiles:
    assert sars2_spike_predictor_phenos_outfile in sars2_spike_predictor_phenos_config
    sars2_spike_predictor_phenos_config[sars2_spike_predictor_phenos_outfile] = os.path.join(
        "results/SARS2-spike-predictor-phenos",
        os.path.basename(sars2_spike_predictor_phenos_config[sars2_spike_predictor_phenos_outfile])
    ) 
sars2_spike_predictor_phenos_infiles = []
for ref_clade_d in sars2_spike_predictor_phenos_config["mutation_phenotype_csvs"].values():
    for phenotype_d in ref_clade_d.values():
        csv_file = os.path.join("SARS2-spike-predictor-phenos", phenotype_d["csv"])
        phenotype_d["csv"] = csv_file
        sars2_spike_predictor_phenos_infiles.append(csv_file)

rule run_SARS2_spike_predictor_phenos:
    """Run the ``SARS2-spike-predictor-phenos`` submodule."""
    input:
        sars2_spike_predictor_phenos_infiles,
        nb="SARS2-spike-predictor-phenos/SARS2-spike-predictor-phenos.ipynb",
    params:
        config_yaml=yaml.round_trip_dump(sars2_spike_predictor_phenos_config),
    output:
        nb="results/SARS2-spike-predictor-phenos/SARS2-spike-predictor-phenos.ipynb",
        config_yaml="results/SARS2-spike-predictor-phenos/config.yaml",
        **{
            outfile: sars2_spike_predictor_phenos_config[outfile]
            for outfile in sars2_spike_predictor_phenos_outfiles
        },
    conda:
        "SARS2-spike-predictor-phenos/environment.yml",
    log:
        log="results/logs/run_SARS2_spike_predictor_phenos.txt",
    shell:
        """
        echo "{params.config_yaml}" > {output.config_yaml} 2> {log}
        papermill -p config_yaml {output.config_yaml} {input.nb} {output.nb} &>> {log}
        """

other_target_files += rules.run_SARS2_spike_predictor_phenos.output


rule compare_natural:
    """Compare DMS (or other) phenotype measurements to natural sequence evolution."""
    input:
        nb="notebooks/compare_natural.ipynb",
        murrell_growth_rates_csv="MultinomialLogisticGrowth/model_fits/rates.csv",
        clade_phenotypes_csv=rules.run_SARS2_spike_predictor_phenos.output.clade_phenotypes_csv,
        clade_phenotypes_randomized_csv=rules.run_SARS2_spike_predictor_phenos.output.clade_phenotypes_randomized_csv,
    output:
        nb="results/notebooks/compare_natural.ipynb",
        pair_growth_dms_csv="results/compare_natural/clade_pair_growth.csv",
        clade_growth_dms_csv="results/compare_natural/clade_growth.csv",
        pair_corr_html="results/compare_natural/clade_pair_growth.html",
        clade_corr_html="results/compare_natural/clade_growth.html",
        pair_ols_html="results/compare_natural/ols_clade_pair_growth.html",
    log:
        log="results/logs/compare_natural.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p murrell_growth_rates_csv {input.murrell_growth_rates_csv} \
            -p clade_phenotypes_csv {input.clade_phenotypes_csv} \
            -p clade_phenotypes_randomized_csv {input.clade_phenotypes_randomized_csv} \
            &> {log}
        """


rule non_rbd_binding_natural:
    """Look at non-RBD mutation effects on ACE2 binding in natural viruses."""
    input:
        dms_summary_csv="results/summaries/summary.csv",
        nb="notebooks/non_rbd_binding_natural.ipynb",
    output:
        nb="results/notebooks/non_rbd_binding_natural.ipynb",
        pango_consensus_seqs_json="results/non_rbd_binding_natural/pango-consensus-sequences_summary.json",
    params:
        pango_consensus_seqs_json="https://raw.githubusercontent.com/corneliusroemer/pango-sequences/71f1c4f249e2a6bc1d2a6b91a310ac1b610cd776/data/pango-consensus-sequences_summary.json",
        yaml=lambda wc, input, output: yaml.round_trip_dump(
            {
                "pango_consensus_seqs_json": output.pango_consensus_seqs_json,
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
        """
        curl {params.pango_consensus_seqs_json} -o {output.pango_consensus_seqs_json} &> {log}
        papermill {input.nb} {output.nb} -y '{params.yaml}' &>> {log}
        """


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
    "DMS measurements versus natural clade growth": {
        "Correlation of change in clade growth versus phenotype for clade pairs":
            rules.compare_natural.output.pair_corr_html,
        "OLS change in clade growth versus phenotype for clade pairs":
            rules.compare_natural.output.pair_ols_html,
        "Correlation of absolute clade growth versus phenotype":
            rules.compare_natural.output.clade_corr_html,
        "Notebook comparing change in clade growth to change in phenotype":
            rules.compare_natural.output.nb,
        "CSV with data for comparison of changes in growth vs phenotype for clade pairs":
            rules.compare_natural.output.pair_growth_dms_csv,
        "CSV with phenotypes and growth for all individual clades":
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
