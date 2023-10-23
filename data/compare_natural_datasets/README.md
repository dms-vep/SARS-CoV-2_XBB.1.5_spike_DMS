# Datasets of effects of mutations for comparison to natural evolution

[EVEscape_XBB_single_mutation_predictions.csv](EVEscape_XBB_single_mutation_predictions.csv) is the data for the single mutation effects to XBB as estimated by EVEscape. This is the file downloaded on Oct-22-2023 from [https://evescape.org/data](https://evescape.org/data) for XBB using the `Mutation Data for Omicron(XBB)` box (which is [this download link](https://api.evescape.org/download_variant_data?curr-virus=COVID19&curr-variant-or-id=Omicron(XBB))). The downloaded file was then edited to remove the `mutations` column (which is very long) and rename the `wt`, `pos`, and `mut` columns to be called `wildtype`, `site`, and `mutant`. Then the various metrics (`EVEscape` etc) were averaged across each site / wildtype / mutant.

[rand_EVEscape.csv](rand_EVEscape.csv) is a version of [EVEscape_XBB_single_mutation_predictions.csv](EVEscape_XBB_single_mutation_predictions.csv) where the `EVEscape` values have been randomized among mutations.

[incremental_Hamming_distance_from_Wuhan-Hu-1.csv](incremental_Hamming_distance_from_Wuhan-Hu-1.csv) just gives a value of one to every mutation that is different from the Wuhan-Hu-1 spike, taken from [QHD43416.1](https://www.ncbi.nlm.nih.gov/protein/QHD43416.1).
