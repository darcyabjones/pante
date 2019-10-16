# PanTE

A pipeline for predicting and masking transposable elements in multiple genomes.

This pipeline is currently mostly developed for personal use, some of the default parameter combinations or tools may be specific for the kinds of organisms I work on (plant-pathogenic Fungi).
I do make an effort to make it usable more generally though.


## What does this do?

PanTE takes a population of genomes and runs several repeat, transposable element, and non-coding RNA prediction tools and merges the results to yield a reasonably comprehensive picture of repeats in your genomes.
My intended use case is for multiple genomes from the same species, but I suppose you could do closely related organisms in the same run too.

To run PanTE you'll just need your genomes and ideally a copy of the RepBase repeat masker formatted database.

The pipeline follows these main steps:

1. Predict non-coding RNA elements using tRNAScan-SE, Infernal (searching against Rfam), and optionally RNAmmer.
2. Predict transposable elements using RepeatModeler, LtrHarvest/LTRDigest, EAHelitron, MiteFinder 2, and MMSeqs2 profile searches against GyDB, selected Pfam models, and a custom set of TE proteins.
3. Combine all TE predictions (except LTRDigest/Harvest) and cluster them to form conservative families using vsearch.
4. Filter the families based on minimum abundance within each genome and presence across the population.
5. Compute multiple sequence alignments for the families using DECIPHER.
6. Classify the families using RepeatClassifier (part of RepeatModeler).
7. Search all genomes for more distant matches to the families (and optionally species models from RepBase/DFAM) using RepeatMasker.
8. Combine all TE and ncRNA predictions into a final GFF file and soft-mask the genomes using this combined set.


The reason that LTRDigest/LTRHarvest is currently excluded from the family clustering is that
the predicted LTRs tend to be quite big or contain nested elements, which tends to group non-related TEs into single clusters.
Eventually I might implement a method to split the LTR predictions up to avoid this issue, but RepeatModeler tends to pick up important parts of LTRs anyway.


## How to not use this pipeline.

There are a couple of pipelines that do repeat annotation, but I haven't seen any that handle multiple genomes particularly well.

Here are some honourable mentions:

- [REPET](https://urgi.versailles.inra.fr/Tools/REPET) is very comprehensive but famously buggy and difficult to install/configure.
- [EDTA](https://github.com/oushujun/EDTA) looks fairly promising and is probably a good choice for Plant genomes.
- [PiRATE](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4763-1) is quite comprehensive.
  It is distributed as a virtual machine, and is run via [Galaxy](https://usegalaxy.org/) within that VM.
  This is probably convenient for people that only have a few genomes to run and would prefer to avoid the command line.


Other pipelines tend to be focussed on inferring repeats from raw reads (e.g. [RepeatExplorer](http://repeatexplorer.org/) and don't offer much for genome annotation.
There is another category of TE pipeline that focusses on insertion site prediction, including [McClintock](https://github.com/bergmanlab/mcclintock), [TEA](http://compbio.med.harvard.edu/Tea/), and [STEAK](https://github.com/applevir/STEAK).
These pipelines are really only useful for organisms with existing, well-curated repeat families and for enabling specific TE-focussed sequencing experiments.


## Quick start

Assuming you have singularity and nextflow installed (See INSTALL).
Say you have a bunch of genome fasta files in a folder `genomes/*.fasta`.
You have downloaded the [RepBase](https://www.girinst.org/repbase/) repeat masker database and the corresponsing RepeatMasker metadata file.

```bash
nextflow run darcyabjones/pante -profile singularity -resume \
  --genomes "genomes/*.fasta" \
  --repbase "downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz" \
  --rm_meta "containers/downloads/RepeatMaskerMetaData-20181026.tar.gz"
```

Will run the full pipeline (except for RNAmmer and RepeatMasker using the "species" model).
Additional databases will be downloaded as part of the pipeline, but you can also download them beforehand and provide them as an argument.

The results will be written to the `results` folder.

## Install

The pipeline is written in [Nextflow](https://www.nextflow.io/), which you will need to install on the executing computer.
The pipeline itself has many dependencies and I have customised RepeatMasker/RepeatModeler a bit, so I HIGHLY recommend that you use the docker or singularity containers.
If you really want to install the software yourself look in the `containers` folder and follow the `Dockerfiles`.

To run the containers, you'll need to install either [Singularity](https://sylabs.io/docs/) (recommended) or [Docker](https://www.docker.com/).

The pipeline itself will pull the containers for you.

On an ubuntu server, the process to install nextflow and singularity might look like this.

```bash
set -eu

sudo apt-get update
sudo apt-get install -y \
    default-jre-headless \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git

VERSION=1.12
OS=linux
ARCH=amd64
wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz
sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz
rm go$VERSION.$OS-$ARCH.tar.gz

echo 'export PATH=/usr/local/go/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

VERSION=3.4.0
wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
tar -xzf singularity-${VERSION}.tar.gz
cd singularity
./mconfig
make -C builddir
sudo make -C builddir install

cd ..
rm -rf -- singularity

curl -s https://get.nextflow.io | bash

./nextflow run darcyabjones/pante --help
```

Because RNAmmer has a restricted license, you'll need to [download the source files yourself](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?rnammer) and build a special container that includes it.
There are instructions for doing this in the [containers/README.md](containers/README.md) file.


## Profiles

The strength of pipeline engines like nextflow is that you can run it on different compute systems
simply by switching some configuration files.

Some preset config files are included in this repo.
You can view these config files in the `conf` directory.

The configuration to use at runtime is controlled by the `-profile` parameter.

Multiple profiles can be specified by separating them with a comma e.g. `-profile laptop,singularity`.
PanTE generally has a separate config file for a compute environment (e.g. cloud, HPC, laptop), and for a software environment (e.g. singularity, docker, local).
It's likely that you'll have to tailor the compute configuration, but you shouldn't need to change the software config so this allows you to mix-and-match.

For more info on configuration see the [nextflow documentation](https://www.nextflow.io/docs/latest/config.html).
You can also raise an issue on the github repository and I'll try to help.


## Parameters

| parameter | default | description |
| :---      | :---    | :---        |
| `--genomes` | Required | A glob of the fasta genomes to search for genes in. The basename of the file is used as the genome name. |
| `--outdir` | results | Where to store the results. |
| `--repbase` | | |
| `--rm_meta` | | |
| `--dfam_hmm` | | |
| `--dfam_hmm_url` | | |
| `--dfam_embl` | | |
| `--dfam_embl_url` | | |
| `--rm_repeatpeps` | | |
| `--rm_species` | | |
| `--mitefinder_profiles` | | |
| `--noinfernal` | | |
| `--rfam` | | |
| `--rfam_url` | | |
| `--rfam_clanin` | | |
| `--rfam_gomapping` | | |
| `--rfam_gomapping_url` | | |
| `--rnammer` | | |
| `--pfam` | | |
| `--gypsydb` | | |
| `--gypsydb_url` | | |
| `--pfam_ids` | `data/pfam_ids.txt` | |
| `--protein_families` | `data/proteins/families.stk` | |
| `--infernal_max_evalue` | 0.00001 | |
| `--mmseqs_max_evalue` | 0.001 | |
| `--min_intra_frequency` | 3 | |
| `--min_inter_proportion` | 0.05 | |
| `--eahelitron_three_prime_fuzzy_level` | | |
| `--eahelitron_upstream_length` | 3000 | |
| `--eahelitron_downstream_length` | 50 | |
| `--ltrharvest_similar` | 85 | |
| `--ltrharvest_vic` | 10 | |
| `--ltrharvest_seed` | 20 | |
| `--ltrharvest_minlenltr` | 100 | |
| `--ltrharvest_maxlenltr` | 7000 | |
| `--ltrharvest_mintsd` | 4 | |
| `--ltrharvest_maxtsd` | 6 | |
| `--mitefinder_threshold` | 0.5 | |
| `--trans_table` | 1 | |

## Exit codes

I'm hoping to add some better error-handling in the future to provide more useful/nextflow-agnostic tips to users.
In the meantime, it's just input parameter validation that is handled elegantly.

- 0: All ok.
- 1: Incomplete parameter inputs.
