# Translation-initiation-associated sequence architecture

## Status
Working reproducibility package. Public deposition has not yet been completed.

## Contents
- scripts/: analysis and figure-generation scripts
- data/orthology/: orthogroup, RBH, and orthology audit tables
- data/upstream_windows/: matched upstream sequences, corrected motif features,
  and null-model outputs
- results/gee_publication/: statistical tables and publication figures
- figures/: copies of figures and associated source tables
- metadata/: genome validation and accession manifest
- environment/: software-version records

## Interpretation
The integrated sequence score and heuristic candidate-motif status are different
analysis endpoints. Neither directly measures translation efficiency, protein
abundance, metabolic flux, or ethanol production.

The dinucleotide-preserving null model preserves the tested sequence constraints,
but uniform sampling over the constrained sequence space has not been established.
Empirical P values from 1,000 replicates have finite simulation resolution.

## Genome accessions
Complete and verify metadata/genome_accession_manifest.tsv against the exact
genome files used in the analysis.

## Reproduction
Review the scripts and document their execution order, dependencies, input files,
and output paths before publication.

Repository DOI/URL: NOT YET ASSIGNED
Public deposition: NOT YET COMPLETED
