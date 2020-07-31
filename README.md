# MAC pipeline : Mapping, Analysis, Comparison
# Installation
## Dependencies
- snakemake (conda install)
- python modules: intervaltree, kerneltree, tqdm, pybedtools, pysam v0.8.4+
- Stringtie
- Gffcompare
- paftools.js : latest version available on the minimap2 github (the conda one isn't working properly).
- bed12ToGtf script
## Installation process
Clone the git repo via git clone or use a release.
# How to use
First edit the `config.yaml` file to set paths to softs, datasets and the name of the cell line used.
Then issue the command:
`snakemake --use-conda -s Snakefile.py --cores [max threads] --cluster "sbatch --cpus-per-task={threads} --mem={resources.ram}" -j 5`
- You can replace "Snakefile.py" with the path to the file but beware that the config file is expected to be in the current working dir.
- The "cores" option set the max threads available to the rules.
- Note that all rules are launched separately with their own RAM and CPU on the cluster : launching snakemake with the defaults values for the cluster is perfectly fine.
