This repository contains a docker image to use CMG-Biotools. 

https://services.healthtech.dtu.dk/services/CMG-biotools/CMGbiotools-2.0/Exercises-CMGbiotools-2.0.pdf

## How to use in huginn

Inside huginn, you can just run the following command:

```bash
conda activate snakemake
singularity shell /software/singularities/cmg_biotools.sif
```

## How to use

First, pull the image (more information in https://github.com/AU-ENVS-Bioinformatics/cmg-biotools/pkgs/container/cmg-biotools)
```bash
docker pull ghcr.io/au-envs-bioinformatics/cmg-biotools:latest
```
Then, run the image (interactive mode):

```bash
docker run -it ghcr.io/au-envs-bioinformatics/cmg-biotools:latest
```

Or bind your current directory, for example: 

```bash
docker run --user $UID:$GID -v $(pwd):/root -it ghcr.io/au-envs-bioinformatics/cmg-biotools:latest
```

## Development instructions

To build the docker image, run the following command:
```bash
docker build -t cmg_biotools .
```
To run the docker image, run the following command:
```bash
docker run -it cmg_biotools
```
