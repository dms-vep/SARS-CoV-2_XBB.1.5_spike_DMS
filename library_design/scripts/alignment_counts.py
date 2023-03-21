"""Get alignment counts."""

import pandas as pd

import requests

table_url = snakemake.params.table_url
alignment_counts = snakemake.output.alignment_counts

print(f"mutation count dataframe from {table_url}")
html_text = requests.get(table_url, verify=False).content.decode("utf-8")
df = pd.read_html(html_text, encoding="utf-8")[0]

print(f"parsing mutation count dataframe ")
cols_of_interest = {
    "WildtypeAA": "wildtype",
    "Position": "site",
    "MutatedAA": "mutant",
    "#Occurrence": "count",
}

valid_mutants = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "-",
    "*",
]

mut_counts = (
    df
    .query("MutatedAA.notnull()", engine="python")
    .query("WildtypeAA.notnull()", engine="python")
    .query("Position.notnull()", engine="python")
    .query("Protein == 'Spike'")
    .rename(columns=cols_of_interest)
    .assign(mutant=lambda x: x["mutant"].map(lambda m: "-" if m == "del" else m))
    .query("mutant in @valid_mutants")
    [cols_of_interest.values()]
    .assign(site=lambda x: x["site"].astype(int))
    .sort_values("count", ascending=False)
    .reset_index(drop=True)
)

print(f"Saving to {alignment_counts=}")
mut_counts.to_csv(alignment_counts, index=False)
