"""Implements ``snakemake`` rule `spatial_distances`"""


import sys

import polyclonal.pdb_utils


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Calculating distances using {snakemake.params.target_chains=}")
spatial_distances = polyclonal.pdb_utils.inter_residue_distances(
    snakemake.input.pdb,
    target_chains=snakemake.params.target_chains,
    target_atom=None,
)

print(f"Writing distances to {snakemake.output.csv=}")
spatial_distances.to_csv(snakemake.output.csv, index=False, float_format="%.5g")