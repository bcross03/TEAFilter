# TEFilter
### Filtering Alignments to Transposable Elements

Scripts for filtering and visualising alignments to transposable elements

**Inputs**
- .fasta file you would like to annotate
- AMRFinderPlus (version >=4.2.5) tabular output for sequence you would like to annotate
- MGE database (or use the one in this repository*)
- Directory of .fasta files for elements to add to database (can be empty but must exist)

**Scripts**
1. Annotation_Filter_Public.py - runs BLAST and filters outputs
2. Annotation_Visualisation_Public.py - Takes outputs from script 1 and generates image of remaining alignments for to help with final manual curation


\* Database is a slightly altered MGEdb, [available here](https://bitbucket.org/mhkj/mgedb/src/master/mgedb/data/sequences.d/mge_records.fna)
