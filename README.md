# PanTE

A pipeline for predicting and masking transposable elements in multiple genomes.

Please note, this is a work in progress.
At the moment it is mainly to support repeat masking for a gene prediction project called [panann](https://github.com/darcyabjones/panann).

The idea is to:

1) take a population of genomes
2) predict TEs/Repeat families in each of them
3) combine and filter the predictions into a non-redundant set
4) transfer consensus annotations onto genomes


There are plenty of tools that can do something like this for individual genomes.
But they tend to be a bit hacky (Perl hacky), parallelism isn't available (or when it is requires horrendous amounts of configuration), often can't take advantage of containers, and tend to be highly specialised for finding plant TEs.


Watch this space!
