"""Implements ``snakemake`` rule to process EVEscape data."""

import sys
import tempfile
import urllib
import zipfile

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

zip_url = snakemake.params.zip_url

# get the spike EVEscape data from the supplementary file ZIP
with tempfile.NamedTemporaryFile() as f:
    print(f"Downloading {zip_url} to temporary file")
    urllib.request.urlretrieve(zip_url, f.name)
    with zipfile.ZipFile(f.name) as zfile:
        with zfile.open(
            "Supp_Table_6_evescape_predictions/spike_evescape_predictions.csv",
            mode="r",
        ) as csvfile:
            df = pd.read_csv(csvfile)
print(f"Read data frame of {len(df)} lines")

print(f"Writing data to {snakemake.output.csv}")

(
    df
    .rename(columns={"i": "site", "wt": "wildtype", "mut": "mutant"})
    .to_csv(snakemake.output.csv, index=False, float_format="%.5g")
)
