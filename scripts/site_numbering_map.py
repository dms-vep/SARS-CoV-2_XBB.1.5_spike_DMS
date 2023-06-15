"""Implements ``snakemake`` rule to map sequential to reference site numbers."""


import string
import subprocess
import sys

import Bio.SeqIO

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# snakemake variables
numbering_reference_accession = snakemake.params.numbering_reference_accession
reference_fasta = snakemake.output.reference
alignment_fasta = snakemake.output.alignment
to_align_fasta = snakemake.output.to_align
prot_fasta = snakemake.input.prot
site_map_csv = snakemake.output.site_numbering_map

print(f"Getting {numbering_reference_accession=} to {reference_fasta=}")
res = subprocess.run(
    [
        "efetch",
        "-format",
        "fasta",
        "-db",
        "protein",
        "-id",
        str(numbering_reference_accession),
    ],
    capture_output=True,
    text=True,
)
with open(reference_fasta, "w") as f:
    f.write(res.stdout)
reference = str(Bio.SeqIO.read(reference_fasta, "fasta").seq)

prot = str(Bio.SeqIO.read(prot_fasta, "fasta").seq)

# align sequences
with open(to_align_fasta, "w") as f:
    f.write(f">reference\n{reference}\n>sequential\n{prot}\n")
res = subprocess.run(
    ["mafft", "--anysymbol", "--amino", to_align_fasta],
    capture_output=True,
    text=True,
    check=True,
)
assert res.returncode == 0
alignment = res.stdout
print(f"Writing alignment to {alignment_fasta}")
with open(alignment_fasta, "w") as f:
    f.write(alignment)

# get the site map
aligned_seqs = list(Bio.SeqIO.parse(alignment_fasta, "fasta"))
assert len(aligned_seqs) == 2
aligned_ref = str(aligned_seqs[0].seq)
aligned_seq = str(aligned_seqs[1].seq)
assert len(aligned_ref) == len(aligned_seq)
records = []
ref_i = seq_i = ref_i_suffix = ref_i_suffix_index = 0
for ref_aa, seq_aa in zip(aligned_ref, aligned_seq):
    if ref_aa != "-":
        ref_i += 1
        ref_i_suffix = ""
        ref_i_suffix_index = 0
    else:
        ref_i_suffix = string.ascii_lowercase[ref_i_suffix_index]
        ref_i_suffix_index += 1
    if seq_aa != "-":
        seq_i += 1
        records.append((seq_i, seq_aa, f"{ref_i}{ref_i_suffix}", ref_aa))
site_map = pd.DataFrame(
    records,
    columns=["sequential_site", "sequential_wt", "reference_site", "reference_wt"],
)
print(f"Writing site map to {site_map_csv}")
site_map.to_csv(site_map_csv, index=False)
