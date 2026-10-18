# Release checklist

## Completed local checks
- [x] Six genome entries recorded in a tab-delimited manifest.
- [x] Local RefSeq annotation names recorded for all six assemblies.
- [x] Standardized genome inputs validated in the analysis workspace.
- [x] Upstream windows mapped to unique target proteins and locus tags.
- [x] All 2,088 upstream windows passed sequence re-extraction validation.
- [x] Corrected motif calls audited against the corrected feature table.
- [x] Matched CCM outputs and the revised five-CCM motif table included.
- [x] Package checksums generated and verified with zero failures.
- [x] Unexpected terminal capture preserved outside the release package.

## Before public deposition
- [x] Confirm all replicon accession versions against the local GenBank records.
- [ ] Review README and reproduction workflow against the final package contents.
- [ ] Test the workflow from a clean environment or document any untested steps.
- [ ] Confirm software versions and commands used for the final analysis.
- [ ] Review the manuscript, figure labels, captions, and reported statistics against final output tables.
- [ ] Remove temporary, redundant, or irrelevant files from the package.
- [ ] Regenerate SHA256SUMS.txt after all final changes and verify every checksum.
- [ ] Create a versioned archive and deposit it in a suitable repository.
- [ ] Add the permanent repository DOI/URL after deposition.

## Interpretation limits
- Motif calls are heuristic sequence-feature classifications, not direct measurements of translation.
- Orthogroup membership is sequence-based and does not by itself establish functional equivalence.
- Dinucleotide-shuffle results should be described with the null-model and sampling limitations stated in the manuscript.
