# pairwise-colocalization-smk

This is a snakemake pipeline for performing colocalization using the coloc R package. It is adapted from an array job pipeline (original scripts in slurm).

The pipeline requires a file containing a list of variants across different summary statistics files and can be passed to the pipeline using ```snakemake --config input_signals="FILEPATH"```.  
