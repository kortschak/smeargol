# Example

This tiny example shows how the GO DAG is painted with a couple of expression data sets.

When run with the following command, the `small.dot` file is created which can then be viewed with GraphViz's `dot` command.

```
smeargol -debug -in sample.counts.gz -map gene_to_go.nt.gz -ontology go.owl.gz > small.dot
```

The graph generated is shown in the README of the smeargol repository.