# AF_MPRA_pipeline

To create the snakemake environment, make sure you're in base conda, then from the `AF_MPRA_pipeline` directory run the following command:

```
conda env create --file snakemake_env.yaml
conda activate snakemake
```

To run the pipeline, with the snakemake conda environment activated, go to the `pipeline` directory and enter the following:

```
snakemake --profile profile
```

You can use the --dry-run option to make sure everything looks correct first without actually running the pipeline.