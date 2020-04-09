# maffer

### extract MSAs from variation graphs

_converts (sorted) genome variation graphs (in GFA or .xg format) to multiple alignment format (MAF)_

## overview

This tool projects between [pangenomic variation graphs](https://pangenome.github.io/) stored in [GFAv1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md), which can be used to encode whole genome alignments, and the multiple alignment format [MAF](http://www.bx.psu.edu/~dcking/man/maf.xhtml), which represents only the linearizable components of such an alignment graph.

Its goal is to allow tools like [seqwish](https://github.com/ekg/seqwish), which efficiently construct whole genome alignment graphs from collections of sequences and pairwise alignments between them, to be applied to comparative genomics problems.

Graphs in GFA must be encoded as is standard in variation graph based methods.
We use the `S`, `L`, and `P` records.
`S` encodes sequences, whose identifiers should be numeric values.
`L` encodes links between sequences.
`P` describes the paths of genomes through the graph.
Graphs of this format are produced by seqwish.

The graph should be sorted in a reasonable way, which can be accomplished with [`odgi sort`](https://pangenome.github.io/odgi/odgi_docs.html#_odgi_sort1).
This is because the blocks in the emitted MAF format are linear components in the sorted graph.
We determine blocks by finding regions in which all path positions versus the sorted graph space are monotonically increasing or decreasing.
In short, structural variation and rearrangement will break the output blocks, as these cannot be represented with a gapped multiple alignment.

## usage

We use cmake to build the tool:

```
cmake -H. -Bbuild && cmake --build build -- -j4
```

Conversion to MAF occurs in two steps.
First, we project the GFA into .xg format:

```
maffer -o graph.xg -g graph.gfa
```

Then, we process the .xg formatted graph (which contains a positional index that we need) to obtain MAF:

```
maffer -i graph.xg >aln.maf
```

## license

MIT

## authors

Erik Garrison
Simon Heumos
