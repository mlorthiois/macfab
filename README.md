# MAC pipeline : Mapping, Analysis, Comparison

# Installation
## Dependencies
    - Snakemake (needs to be installed through conda)
    - (All tools are installed and maintained by Snakemake)

## Installation process
Clone the git repo via `git clone git@gitlab.com:bioinfog/macfab.git` (or `git clone https://gitlab.com/bioinfog/macfab.git`).

# How to use
1. Edit the `config.yaml` file to set paths to datasets and metadata.

2. Enter in the repo, load snakemake env and issue the command:
```
. /local/env/envsnakemake-5.20.1.sh

snakemake --use-conda --cores [number] -j 5
```

# Informations
- The "cores" option set the max threads available to the rules.
- Your gtf, fa and fatsq (reported in `config.yaml`) can be gzipped or not.
- You can add other parameters in the command line ([description here](https://snakemake.readthedocs.io/en/stable/executing/cli.html)).
- You visualize via a [DAG](https://en.wikipedia.org/wiki/Directed_acyclic_graph) using `snakemake --dag -n ` and then [dedicated website](https://dreampuf.github.io/GraphvizOnline) to visualize it.
