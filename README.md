# TEAFilter
### Transposable Element Alignment Filter

Scripts for filtering and visualising alignments to transposable elements

**Inputs**
- .fasta file you would like to annotate
- AMRFinderPlus (version >=4.2.5) tabular output
- Multi-fasta MGE database*
- Directory of .fasta files for elements to add to database (can be empty but must exist, see Additional_Elements/ in this repository for example)
- IntegronFinder output (optional)

**Scripts**
1. Annotation_Filter.py - runs BLAST and filters outputs
2. Annotation_Visualisation.py - Takes outputs from script 1 and generates image of remaining alignments to help with final manual curation


\* MGEdb contains many transposons, insertion sequences and ICE, [available here](https://bitbucket.org/mhkj/mgedb/src/master/mgedb/data/sequences.d/mge_records.fna) \
\* TnCentral contains well-curated transposons with some insertion sequences and integrons, [available here](https://tncentral.ncc.unesp.br/data_download)
