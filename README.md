# smeargol

`smeargol` is a tool for non-redundantly assigning gene count data to Gene Ontology terms associated with the genes. It is based on ideas from [Fruzangohar _et al._](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0170486).

`smeargol` distributes count data across the Gene Ontology DAG provided and prints the GO terms, their roots and depths and distributed counts in a tsv table to stdout. It logs gene identifiers that do not have GO term annotations to stderr. The graph analysis assumes Ensembl gene identifiers and Gene Ontology graph structure.

The input counts file is a tab-delimited file with the first column being Ensembl gene ID (ENSG00000000000) and remaining columns being count data. The first row is expected to be labelled with the first column being Geneid and the remaining columns holding the names of the samples.

The Gene Ontology is required to be in Owl format. The file can be obtained from http://current.geneontology.org/ontology/go.owl.

The ENSG to GO mapping is expected to be in RDF N-Triples or N-Quads in the form:

```
<obo:GO_0000000> <local:annotates> <ensembl:ENSG00000000000> .
```

for each GO term to Ensembl gene annotation.

All input files are expected to be gzip compressed and the output is written uncompressed to standard output and standard error.