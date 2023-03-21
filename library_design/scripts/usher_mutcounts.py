"""Get mutation counts from translated UShER tree."""

import pandas as pd


translated_muts_tsv = snakemake.input.translated_muts_tsv
mut_counts_csv = snakemake.output.mut_counts_csv

df = (
    pd.read_csv(translated_muts_tsv, sep="\t")
    .assign(mutation=lambda x: x["aa_mutations"].str.split(";"))
    [["mutation"]]
    .explode("mutation")
    .query("mutation.str.startswith('S:')")
    .assign(
        mutation=lambda x: x["mutation"].str[2: ],
        site=lambda x: x["mutation"].str[1: -1].astype(int),
        mutant_aa=lambda x: x["mutation"].str[-1],
    )
    .groupby(["site", "mutant_aa"], as_index=False)
    .aggregate(count=pd.NamedAgg("site", "count"))
    .sort_values("count", ascending=False)
    .to_csv(mut_counts_csv, index=False)
)
