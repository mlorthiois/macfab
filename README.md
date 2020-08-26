# MAC pipeline : Mapping, Analysis, Comparison
# Installation
## Dependencies
- snakemake (needs to be installed through conda)
- python modules: intervaltree, kerneltree, tqdm, pybedtools, pysam v0.8.4+
- talon (**pip install** or locate softs in `~/.local/bin/`)
- Stringtie : https://github.com/gpertea/stringtie
- Gffcompare : http://ccb.jhu.edu/software/stringtie/gffcompare.shtml
- paftools.js : latest version available on the minimap2 github (the conda one isn't working properly).
- bed12ToGtf script
## Installation process
Clone the git repo via git clone or use a release.
# How to use
First edit the `config.yaml` file to set paths to softs, datasets and metadata.
Then issue the command:
`snakemake --use-conda -s Snakefile.py --cores [number] --cluster "sbatch --cpus-per-task={threads} --mem={resources.ram}" -j 5`
- You can replace "Snakefile.py" with the path to the file but beware that the config file is expected to be in the current working dir.
- The "cores" option set the max threads available to the rules.
- Note that all rules are launched separately with their own RAM and CPU on the cluster : launching snakemake with the defaults values for the cluster is fine.
- full options description : https://snakemake.readthedocs.io/en/stable/executing/cli.html
