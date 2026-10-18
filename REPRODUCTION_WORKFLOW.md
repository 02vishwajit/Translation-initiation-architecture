# Reproduction workflow

## Scope
Comparative analysis of upstream translation-initiation-associated sequence
features in six bacterial genomes. Motif presence is a sequence-based proxy,
not a measurement of translation efficiency, protein abundance, or metabolic flux.

## Included analysis assets
- Standardized genome metadata and validation report
- Orthology and single-copy orthogroup audit tables
- Extracted 50-nt upstream windows and corrected sequence-feature tables
- Dinucleotide-preserving null-model outputs and validation reports
- Adjusted GEE results, sensitivity tables, and publication figures
- CCM-focused motif refinement and matched-comparison outputs
- Python analysis scripts and software-version record

## Reproduction order
1. Verify the genome accession and annotation manifest against the exact
   genome files used in the analysis.
2. Verify software versions and install the required dependencies.
3. Validate the standardized genome files and protein identifiers.
4. Reproduce orthology inference and the single-copy orthogroup selection.
5. Re-extract and validate upstream windows against the annotated genomes.
6. Run the corrected motif scanner and audit its calls.
7. Run the composition-preserving null models and their validation checks.
8. Refit adjusted GEE models and regenerate tables and figures.
9. Reproduce the CCM-focused comparison separately.
10. Compare regenerated outputs with the supplied result files.

Scripts may contain paths specific to the original workstation. Inspect and
adapt paths before running them on another machine. Do not assume that the
entire workflow is reproducible until it has been executed from a clean
environment and the regenerated outputs have been compared.

## Interpretation limitations
- Motif calls are heuristic sequence classifications.
- A candidate ribosome-binding-site-like sequence does not establish
  translation initiation experimentally.
- Orthogroup membership is a computational inference.
- Genome annotation and strain differences can influence comparisons.
- Empirical null-model P values are limited by the number of replicates.
- No causal claim about ethanol production or translation efficiency follows
  from these sequence analyses alone.

## Release status
This is a working reproducibility package. A public repository URL or DOI
has not yet been assigned. Genome metadata, clean-environment reproduction,
and release checks must be completed before describing it as fully validated.
