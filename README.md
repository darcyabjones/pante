# PanTE2 - version update!!!

An updated version of PanTE. Please run the following command using singularity (nexflow version has to be <=22.10.7)
```
NXF_ANSI_LOG=false nextflow run KristinaGagalova/pante2 -r version-update \
  -profile singularity,standard \
  -resume \
  --genomes "test/*.fasta" \
  --outdir "test/results"
```

Otherwise clone the repo locally and run the command **inside the ```pante2``` directory**
```
cd pante2
NXF_ANSI_LOG=false nextflow run ./main.nf \
  -profile singularity,standard \
  -resume \
  --genomes "test/*.fasta" \
  --outdir "test/results"
```

***
# PanTE

A pipeline for predicting and masking transposable elements in multiple genomes.

This pipeline is currently mostly developed for personal use, some of the default parameter combinations or tools may be specific for the kinds of organisms I work on (plant-pathogenic Fungi).
I do make an effort to make it usable more generally though.


## What does this do?

PanTE takes a population of genomes and runs several repeat, transposable element, and non-coding RNA prediction tools and merges the results to yield a reasonably comprehensive picture of repeats in your genomes.
My intended use case is for multiple genomes from the same species, but I suppose you could do closely related organisms in the same run too.

To run PanTE you'll just need your genomes and optionally a copy of the [RepBase repeat masker formatted database](https://www.girinst.org/server/RepBase/index.php).

The pipeline follows these main steps:

1. Predict non-coding RNA elements using [tRNAScan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/), [Infernal](http://eddylab.org/infernal/) (searching against [Rfam](https://rfam.xfam.org/)), and optionally [RNAmmer](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=rnammer&version=1.2&packageversion=1.2&platform=Unix).
2. Predict transposable elements using [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/), [LTR](http://genometools.org/tools/gt_ltrharvest.html)[Harvest](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-18)/[LTR](http://genometools.org/tools/gt_ltrdigest.html)[Digest](https://academic.oup.com/nar/article/37/21/7002/1420683), [EAHelitron](https://github.com/dontkme/EAHelitron), [MiteFinder 2](https://github.com/screamer/miteFinder), and [MMSeqs2](https://github.com/soedinglab/MMseqs2) profile searches against [GyDB](http://www.gydb.org/index.php/Main_Page), selected [Pfam](http://www.gydb.org/index.php/Main_Page) models, and a custom set of TE proteins derived from the [TransposonPSI](http://transposonpsi.sourceforge.net/) and [LTR_retriever](https://github.com/oushujun/LTR_retriever/tree/master/database) libraries.
3. Combine all TE predictions (except LTRDigest/Harvest) and cluster them to form conservative families using [vsearch](https://github.com/torognes/vsearch).
4. Filter the families based on minimum abundance within each genome and presence across the population.
5. Compute multiple sequence alignments for the families using [DECIPHER](http://www2.decipher.codes/).
6. Classify the families using RepeatClassifier (part of RepeatModeler).
7. Search all genomes for more distant matches to the families (and optionally species models from [RepBase](https://www.girinst.org/repbase/)/[DFAM](https://dfam.org/home)) using [RepeatMasker](http://www.repeatmasker.org/RMDownload.html).
8. Combine all TE and ncRNA predictions into a final GFF files and soft-mask the genomes using this combined set.


The reason that LTRDigest/LTRHarvest is currently excluded from the family clustering is that
the predicted LTRs tend to be quite big or contain nested elements, which tends to group non-related TEs into single clusters.
Eventually I might implement a method to split the LTR predictions up to avoid this issue, but RepeatModeler tends to pick up important parts of LTRs anyway.


## How to not use this pipeline.

There are a couple of pipelines that do repeat annotation, but I haven't seen any that handle multiple genomes particularly well.

Here are some honourable mentions:

- [REPET](https://urgi.versailles.inra.fr/Tools/REPET) is very comprehensive but famously buggy and difficult to install/configure.
- [EDTA](https://github.com/oushujun/EDTA) looks fairly promising and is probably a good choice for plant genomes.
- [PiRATE](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4763-1) is quite comprehensive.
  It is distributed as a virtual machine, and is run via [Galaxy](https://usegalaxy.org/) within that VM.
  This is probably convenient for people that only have a few genomes to run and would prefer to avoid the command line.


Other pipelines tend to be focussed on inferring repeats from raw reads (e.g. [RepeatExplorer](http://repeatexplorer.org/) and don't offer much for genome annotation.
There is another category of TE pipeline that focusses on insertion site prediction, including [McClintock](https://github.com/bergmanlab/mcclintock), [TEA](http://compbio.med.harvard.edu/Tea/), and [STEAK](https://github.com/applevir/STEAK).
These pipelines are really only useful for organisms with existing, well-curated repeat families and for enabling specific TE-focussed sequencing experiments.


## Running PanTE

Assuming you have singularity and nextflow installed (See INSTALL).
Say you have a bunch of genome fasta files in a folder `genomes/*.fasta`.

```bash
nextflow run KristinaGagalova/pante2 -profile singularity -resume --genomes "genomes/*.fasta"
```

Will run the full pipeline (except for RNAmmer and RepeatMasker using the "species" model).
Additional databases like [Dfam](https://dfam.org/home), [Rfam](https://rfam.xfam.org/), [GyDB](http://www.gydb.org/index.php/Main_Page), and selected [Pfam](https://pfam.xfam.org/) models will be downloaded as part of the pipeline, but you can also download them beforehand and provide them as an argument.

The results will be written to the `results` folder.


If you provide the `--species` parameter, a separate pass of RepeatMasker will be run using the Dfam and/or RepBase databases instead of the custom libraries.
If you would like to include this extra step I would highly recommend providing the [RepBase](https://www.girinst.org/repbase/) repeat masker database and the corresponsing [RepeatMasker metadata files](http://www.repeatmasker.org/libraries/) if you have access.
If you do have access to RepBase, it's probably worth using it even if you aren't using the `--species` option because it might help improve family annotation.

The value given to `--species` can be any NCBI taxonomy name and is provided to the RepeatMasker option `-species`.

```bash
nextflow run KristinaGagalova/pante2 -profile singularity -resume \
  --genomes "genomes/*.fasta" \
  --repbase "downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz" \
  --rm_meta "downloads/RepeatMaskerMetaData-20181026.tar.gz" \
  --species "fungi"
```


If you would like to include [RNAmmer](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=rnammer&version=1.2&packageversion=1.2&platform=Unix) rDNA predictions, you'll need to either install it on all machines that you're running the pipeline on, or you can build the extra container that does the install for you (See [containers/README.md](containers#proprietary-software)).

Then you can provide the `--rnammer` flag to enable those steps.
Here i'm assuming that you've installed RNAmmer locally.

```bash
nextflow run /path/to/pante2/main.nf -profile singularity -resume \
  --genomes "genomes/*.fasta" \
  --species "fungi" \
  --rnammer \
  -with-singularity "containers/singularity/pante2-rnammer.sif"
```
or
```bash
nextflow run /path/to/pante2/main.nf -profile docker -resume \
  --genomes "genomes/*.fasta" \
  --species "fungi" \
  --rnammer \
  -with-docker kristinagagalova/pante2-rnammer:v1.0.0
```


## Install

The pipeline is written in [Nextflow](https://www.nextflow.io/), which you will need to install on the executing computer.
The pipeline itself has many dependencies and I have customised RepeatMasker/RepeatModeler a bit, so I HIGHLY recommend that you use the docker or singularity containers.
If you really want to install the software yourself look in the `containers` folder and follow the `Dockerfiles`.

To run the containers, you'll need to install either [Singularity](https://sylabs.io/docs/) (recommended) or [Docker](https://www.docker.com/).

The pipeline itself will pull the containers for you from [Sylabs Cloud](https://cloud.sylabs.io/library/kristinagagalova/default/pante2) or [DockerHub](https://hub.docker.com/r/kristinagagalova/pante2).

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

./nextflow run KristinaGagalova/pante2 --help
```

Because RNAmmer has a restricted license, you'll need to [download the source files yourself](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=rnammer&version=1.2&packageversion=1.2&platform=Unix) and either install it locally or build a special container that includes it.
There are instructions for doing this [here](containers#proprietary-software).


## Profiles

The strength of pipeline engines like nextflow is that you can run it on different compute systems
simply by switching some configuration files.

Some preset config files are included in this repo.
You can view these config files in the `conf` directory.

The configuration to use at runtime is controlled by the `-profile` parameter.

Multiple profiles can be specified by separating them with a comma e.g. `-profile standard,singularity`.
PanTE generally has a separate config file for a compute environment (e.g. cloud, HPC, laptop), and for a software environment (e.g. singularity, docker, local).
It's likely that you'll have to tailor the compute configuration, but you shouldn't need to change the software config so this allows you to mix-and-match.

Available profiles for containerised software environments are:

- `singularity` - Use a pre-built singularity image containing all non-proprietary software available from https://cloud.sylabs.io/library/kristinagagalova/default/pante2.
- `docker` - Use a pre-build docker container. Like `singularity`. Available from https://hub.docker.com/r/kristinagagalova/pante2.

If you don't specify a software environment profile, it is assumed that all dependencies are installed locally and available on your `PATH`.

NOTE: the docker profiles assume that running `docker` does not require `sudo`.
To use this you'll need to configure docker for "sudo-less" operation ([instructions here](https://docs.docker.com/install/linux/linux-postinstall/)).
If you don't like this, try singularity :)


Available compute profiles are:

- `standard` - (Default) Appropriate for running on a laptop with 4 CPUs and ~8GB RAM.
- `nimbus` - Appropriate for cloud VMs or a local desktop with 2-16 CPUs and ~4-64 GB RAM each (depending on the Nimbus VM flavour you have access to).
- `pawsey_zeus` - Is a config for running on the [Pawsey Zeus](https://pawsey.org.au/systems/zeus/) compute cluster using SLURM.
  Use this as more of a template for setting up your own profile as HPC configuration is pretty specific (and in this case contains some hard coded user options, sorry).


To add your own profile, you can use the files in `./conf` as a template, and make sure you add them to `nextflow.config` under the `profiles` block.

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
| `--rfam` | Optional | Pre-downloaded Rfam CM models (un-gzipped) to use. |
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
| `--min_intra_frequency` | 4 | The minimum number of copies a clustered repeat family must have within a genome for it to be considered "present". |
| `--min_inter_proportion` | 0.2 | The minimum proportion of genomes that the clustered repeat family must be present in (after `--min_intra_frequency`) to be considered a geniune family. |
| `--repeatmodeler_min_len` | 10000 | The minimum scaffold length to allow for predicting repeats in repeatmodeler. Scaffolds smaller than this will be removed to avoid sampling bias. |
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


## Examples

A test dataset and example command is provided in the `test` folder.
On a laptop this takes about an hour to run.


## Exit codes

I'm hoping to add some better error-handling in the future to provide more useful/nextflow-agnostic tips to users.
In the meantime, it's just input parameter validation that is handled elegantly.

- 0: All ok.
- 1: Incomplete parameter inputs.
