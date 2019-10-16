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
| `--outdir` | `results`| The directory to store the results in. |
| `--repbase` | Optional | The RepBase RepeatMasker edition tarball to use to construct the repeatmasker database. Download from [https://www.girinst.org/server/RepBase/index.php](https://www.girinst.org/server/RepBase/index.php). |
| `--rm_meta` | Optional | The RepeatMasker meta tarball to use to construct the repeatmasker database. Download from [http://www.repeatmasker.org/libraries/](http://www.repeatmasker.org/libraries/). Make sure the version matches the version of Repbase if you're using RepBase. |
| `--dfam_hmm` | Optional | Pre downloaded Dfam HMMs to use. Will download latest if this isn't provided. |
| `--dfam_hmm_url` | URL to `Dfam.hmm.gz` | The url to download the Dfam HMMs from if `--dfam_hmm` isn't provided. |
| `--dfam_embl` | Optional | Pre downloaded Dfam consensus sequences to use. Will download latest if this isn't provided. |
| `--dfam_embl_url` | URL to `Dfam.embl.gz` | The url to download the Dfam consensus sequences from if `--dfam_embl` isn't provided. |
| `--rm_repeatpeps` | Optional | Repeat proteins to use for repeatmasker. By default this is taken from the RepeatMasker `Library/RepeatPeps.lib` and assumes that you're using the containers. |
| `--rm_species` | Optional | An NCBI taxonomy name to use to predict transposable elements from RepBase with. Something like `fungi` usually works fine. |
| `--mitefinder_profiles` | Optional | A text file for MiteFinderII containing profiles to search for. Corresponds to [https://github.com/screamer/miteFinder/blob/master/profile/pattern_scoring.txt]. By default will use a file pointed to by the `MITEFINDER_PROFILE` environment variable, which is set in the provided containers. |
| `--noinfernal` | false | Don't run Infernal `cmscan` against Rfam. This can save some time. |
| `--rfam` | Optional | Pre-downloaded Rfam CM models (gzipped) to use. |
| `--rfam_url` | URL to `Rfam.cm.gz` | The url to download Rfam CM models from if `--rfam` isn't provided. Will not download if `--noinfernal` |
| `--rfam_clanin` | Optional | Pre-downloaded Rfam clan information to use. |
| `--rfam_clanin_url` | URL to `Rfam.clanin` | The URL to download Rfam clan info from if `--rfam_clanin` isn't provided. |
| `--rfam_gomapping` | Optional | Pre-downloaded Rfam GO term mappings to use. |
| `--rfam_gomapping_url` | URL to `rfam2go` | The URL to download Rfam GO term mappings from if `--rfam_gomapping` isn't provided. |
| `--rnammer` | false | Run RNAmmer analyses on the genomes. Assumes that you are using the containers with RNAmmer installed or have otherwise set RNAmmer. Will fail if it isn't installed. |
| `--pfam` | Optional | A glob of Pfam stockholm formatted alignments (not gzipped) to use to search against the genomes. |
| `--pfam_ids` | `data/pfam_ids.txt` | A file containing a list of Pfam accessions to download and use if `--pfam` isn't provided. |
| `--gypsydb` | Optional | A glob of stockholm formatted alignments from GyDB to search against the genomes. |
| `--gypsydb_url` | URL to `GyDB_collection.zip` | The URL to download GyDB from if `--gypsydb` is not provided. |
| `--protein_families` | `data/proteins/families.stk` | A stockholm formatted file of custom aligned protein families to search against the genomes. |
| `--infernal_max_evalue` | 0.00001 | The maximum e-value to use to consider `cmscan` matches significant. |
| `--mmseqs_max_evalue` | 0.001 | The maximum e-value to use to consider `mmseqs` profile matches against the genomes significant. |
| `--min_intra_frequency` | 3 | The minimum number of copies a clustered repeat family must have within a genome for it to be considered "present". |
| `--min_inter_proportion` | 0.05 | The minimum proportion of genomes that the clustered repeat family must be present in (after `--min_intra_frequency`) to be considered a geniune family. |
| `--eahelitron_three_prime_fuzzy_level` | 3 | Passed on to the EAHelitron parameter `-r`. |
| `--eahelitron_upstream_length` | 3000 | Passed on to the EAHElitron parameter `-u`. |
| `--eahelitron_downstream_length` | 50 | Passed on to the EAHelitron parameter `d`. |
| `--ltrharvest_similar` | 85 | Passed on to the LTRHarvest parameter `-similar`.  |
| `--ltrharvest_vic` | 10 | Passed on to the LTRHarvest parameter `-vic`. |
| `--ltrharvest_seed` | 20 | Passed on to the LTRHarvest parameter `-seed`. |
| `--ltrharvest_minlenltr` | 100 | Passed on to the LTRHarvest parameter `-minlenltr`. |
| `--ltrharvest_maxlenltr` | 7000 | Passed on to the LTRHarvest parameter `-maxlenltr`. |
| `--ltrharvest_mintsd` | 4 | Passed on to the LTRHarvest parameter `-mintsd`. |
| `--ltrharvest_maxtsd` | 6 | Passed on to the LTRHarvest parameter `-maxtsd`. |
| `--mitefinder_threshold` | 0.5 | Passed on to the MiteFinderII parameter `-threshold`. |
| `--trans_table` | 1 | The ncbi translation table number to use for MMSeqs searches. |


## Exit codes

I'm hoping to add some better error-handling in the future to provide more useful/nextflow-agnostic tips to users.
In the meantime, it's just input parameter validation that is handled elegantly.

- 0: All ok.
- 1: Incomplete parameter inputs.
