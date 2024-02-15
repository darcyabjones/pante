#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # PanTE

    A pipeline for predicting and masking transposable elements in multiple genomes.

    ## Examples

    Say you have a bunch of genomes in `genomes`, you can predict
    ncRNAs and TEs like so:

    ```bash
    nextflow run KristinaGagalova/pante2 \
      -profile singularity \
      -resume \
      --genomes "genomes/*.fasta"
    ```

    If you have access to RepBase you can include that:

    ```
    nextflow run KristinaGagalova/pante2 -profile singularity -resume \
      --genomes "genomes/*.fasta" \
      --repbase "downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz" \
      --rm_meta "downloads/RepeatMaskerMetaData-20181026.tar.gz" \
      --rm_species "fungi"
    ```

    ## Parameters

    --genomes <glob>
        Required
        A glob of the fasta genomes to search for genes in.
        The basename of the file is used as the genome name.

    --outdir <path>
        Default: `results`
        The directory to store the results in.

    --repbase <path>
        Optional
        The RepBase RepeatMasker edition tarball to use to construct
        the repeatmasker database. Download from
        https://www.girinst.org/server/RepBase/index.php.

    --rm_meta <path>
        Optional
        The RepeatMasker meta tarball to use to construct the
        repeatmasker database. Download from
        http://www.repeatmasker.org/libraries/. Make sure the version
        matches the version of Repbase if you're using RepBase.

    --dfam_h5 <glob>
        Optional
        Pre downloaded Dfam h5 to use. Will download latest
        if this isn't provided.

    --dfam_h5_url <url>
        https://www.dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.0.h5.gz
        The url to download the Dfam HMMs from if `--dfam_h5` isn't provided.

    --dfam_embl <path>
        Optional
        Pre downloaded Dfam consensus sequences to use. Will download
        latest if this isn't provided.

    --dfam_embl_url <url>
        http://dfam.org/releases/current/families/Dfam.embl.gz
        The url to download the Dfam consensus sequences from
        if `--dfam_embl` isn't provided.

    --rm_repeatpeps <path>
        Optional
        Repeat proteins to use for repeatmasker. By default this
        is taken from the RepeatMasker `Library/RepeatPeps.lib`
        and assumes that you're using the containers.

    --rm_species <str>
        Optional
        An NCBI taxonomy name to use to predict transposable
        elements from RepBase with. Something like `fungi` usually works fine. |

    --mitefinder_profiles <path>
        Optional
        A text file for MiteFinderII containing profiles to search for.
        Corresponds to https://github.com/screamer/miteFinder/blob/master/profile/pattern_scoring.txt.
        By default will use a file pointed to by the
        `MITEFINDER_PROFILE` environment variable, which is
        set in the provided containers.

    --noinfernal
        false
        Don't run Infernal `cmscan` against Rfam. This can save some time.

    --rfam <path>
        Optional
        Pre-downloaded Rfam CM models (un-gzipped) to use.

    --rfam_url <url>
        ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
        The url to download Rfam CM models from if `--rfam`
        isn't provided. Will not download if `--noinfernal`.

    --rfam_clanin
        Optional
        Pre-downloaded Rfam clan information to use.

    --rfam_clanin_url <url>
        ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
        The URL to download Rfam clan info from if `--rfam_clanin`
        isn't provided.

    --rfam_gomapping <path>
        Optional
        Pre-downloaded Rfam GO term mappings to use.

    --rfam_gomapping_url <url>
        ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/rfam2go/rfam2go
        The URL to download Rfam GO term mappings from if
        `--rfam_gomapping` isn't provided.

    --rnammer
        false
        Run RNAmmer analyses on the genomes. Assumes that you
        are using the containers with RNAmmer installed or
        have otherwise set RNAmmer. Will fail if it isn't installed.

    --pfam <glob>
        Optional
        A glob of Pfam stockholm formatted alignments (not gzipped)
        to use to search against the genomes.

    --pfam_ids <path>
        `data/pfam_ids.txt`
        A file containing a list of Pfam accessions to download
        and use if `--pfam` isn't provided.

    --gypsydb <path>
        Optional
        A glob of stockholm formatted alignments from
        GyDB to search against the genomes.

    --gypsydb_url <url>
        https://gydb.org/extensions/Collection/collection/db/GyDB_collection.zip
        The URL to download GyDB from if `--gypsydb` is not provided.

    --protein_families <path>
        `data/proteins/families.stk`
        A stockholm formatted file of custom aligned protein families
        to search against the genomes.

    --infernal_max_evalue <float>
        0.00001
        The maximum e-value to use to consider `cmscan` matches significant.

    --mmseqs_max_evalue
        0.001
        The maximum e-value to use to consider `mmseqs` profile
        matches against the genomes significant.

    --min_intra_frequency <int>
        4
        The minimum number of copies a clustered repeat family
        must have within a genome for it to be considered "present".

    --min_inter_proportion <float>
        0.2
        The minimum proportion of genomes that the clustered repeat
        family must be present in (after `--min_intra_frequency`)
        to be considered a geniune family.

    --repeatmodeler_min_len <int>
        10000
        The minimum scaffold length to allow for predicting repeats in repeatmodeler.
        Scaffolds smaller than this will be removed to avoid sampling bias.

    --eahelitron_three_prime_fuzzy_level <int>
        3
        Passed on to the EAHelitron parameter `-r`.

    --eahelitron_upstream_length <url>
        3000
        Passed on to the EAHElitron parameter `-u`.

    --eahelitron_downstream_length <int>
        50
        Passed on to the EAHelitron parameter `d`.

    --ltrharvest_similar <int>
        85
        Passed on to the LTRHarvest parameter `-similar`.

    --ltrharvest_vic <int>
        10
        Passed on to the LTRHarvest parameter `-vic`.

    --ltrharvest_seed <int>
        20
        Passed on to the LTRHarvest parameter `-seed`.

    --ltrharvest_minlenltr <int>
        100
        Passed on to the LTRHarvest parameter `-minlenltr`.

    --ltrharvest_maxlenltr <int>
        7000
        Passed on to the LTRHarvest parameter `-maxlenltr`.

    --ltrharvest_mintsd <int>
        4
        Passed on to the LTRHarvest parameter `-mintsd`.

    --ltrharvest_maxtsd <int>
        6
        Passed on to the LTRHarvest parameter `-maxtsd`.
    --mitefinder_threshold <float>
        0.5
        Passed on to the MiteFinderII parameter `-threshold`.

    --trans_table <int>
        1
        The ncbi translation table number to use for MMSeqs searches.

    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}


params.genomes = false

// RepeatMasker config
params.repbase = false
params.rm_meta = false

params.dfam_h5 = false
params.dfam_h5_url = "https://www.dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.0.h5.gz"

params.dfam_embl = false
params.dfam_embl_url = "http://dfam.org/releases/current/families/Dfam.embl.gz"

params.rm_repeatpeps = false
params.rm_species = false

// Mitefinder comes with a set of profiles to use for searching.
// We can fetch this from the installation directory in the containers.
// But if you aren't using the containers or want to use a different
// set of profiles use this.
params.mitefinder_profiles = false

// Infernal takes quite a long time to run and tRNAscan and RNAmmer cover
// a lot of what we want it for. So option is offered to disable it.
params.noinfernal = false
params.rfam = false
params.rfam_url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
params.rfam_clanin = false
params.rfam_clanin_url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"
params.rfam_gomapping = false
params.rfam_gomapping_url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/rfam2go"
params.rnammer = false


params.pfam = false
params.gypsydb = false
params.gypsydb_url = "https://gydb.org/extensions/Collection/collection/db/GyDB_collection.zip"
params.pfam_ids = "$baseDir/data/pfam_ids.txt"
params.protein_families = "$baseDir/data/proteins/families.stk"

params.repeatmodeler_min_len = 1000

params.infernal_max_evalue = 0.00001
params.mmseqs_max_evalue = 0.001

params.min_intra_frequency = 4
params.min_inter_proportion = 0.2

params.eahelitron_three_prime_fuzzy_level = 3
params.eahelitron_upstream_length = 3000
params.eahelitron_downstream_length = 50

params.ltrharvest_similar = 85
params.ltrharvest_vic = 10
params.ltrharvest_seed = 20
params.ltrharvest_minlenltr = 100
params.ltrharvest_maxlenltr = 7000
params.ltrharvest_mintsd = 4
params.ltrharvest_maxtsd = 6

params.mitefinder_threshold = 0.5

params.trans_table = 1

def run_infernal = !params.noinfernal


def exclude_unclean = { filename -> filename.endsWith(".unclean") ? null : filename }


if ( params.genomes ) {
    genomes = Channel
        .fromPath(params.genomes, checkIfExists: true, type: "file")
        .map { file -> [file.baseName, file] }
} else {
    log.info "Hey I need some genomes to assess please."
    exit 1
}


if ( params.repbase ) {
    Channel
        .fromPath(params.repbase, checkIfExists: true, type: "file")
        .first()
        .set { repbase }
} else {
    Channel
        .value( file("REPBASE_WAS_NULL") )
        .set { repbase }
}

if ( params.rm_meta ) {
    Channel
        .fromPath(params.rm_meta, checkIfExists: true, type: "file")
        .first()
        .set { rmMeta }
} else {

    Channel
        .value( file("RM_META_WAS_NULL") )
        .set { rmMeta }
}


if ( params.dfam_h5 ) {

    Channel
        .fromPath(params.dfam_h5, checkIfExists: true, type: "file")
        .first()
        .set { dfamH5 }

} else {

    process getDfamH5 {

        label "download"
        label "small_task"
        time "24h"

        publishDir "${params.outdir}/downloads"

        output:
        file "dfam38.0.h5.gz" into dfamH5

        script:
        """
        wget \
          --no-check-certificate \
	  --tries=1 \
          -c \
          -O dfam38.0.h5.gz \
          "${params.dfam_h5_url}"
        """
    }
}


if ( params.dfam_embl ) {

    Channel
        .fromPath(params.dfam_embl, checkIfExists: true, type: "file")
        .first()
        .set { dfamEMBLs }

} else {

    process getDfamEMBLs {

        label "download"
        label "small_task"
        time "12h"

        publishDir "${params.outdir}/downloads"

        output:
        file "dfam.embl.gz" into dfamEMBLs

        script:
        """
        wget \
          --no-check-certificate \
          -c \
          -O dfam.embl.gz \
          "${params.dfam_embl_url}"
        """
    }
}


if ( params.rm_repeatpeps ) {

    Channel
        .fromPath(params.rm_repeatpeps, checkIfExists: true, type: "file")
        .first()
        .set { rmRepeatPeps }

} else {

    process getRMRepeatPeps {

        label "repeatmasker"
        label "small_task"
        time "20m"

        publishDir "${params.outdir}/downloads"

        output:
        file "RepeatPeps.lib" into rmRepeatPeps

        script:
        """
        cp -L "\${RMASK_PREFIX}/Libraries/RepeatPeps.lib" RepeatPeps.lib
        """
    }
}


if ( params.pfam_ids ) {
    Channel
        .fromPath(params.pfam_ids, checkIfExists: true, type: "file")
        .first()
        .set { pfamIds }
} else {
    pfamLtrIds = Channel.empty()
}


if ( params.rfam && params.rfam_clanin ) {
    Channel
        .value([
            file(params.rfam, checkIfExists: true),
            file(params.rfam_clanin, checkIfExists: true),
        ])
        .set { rfam }

} else if ( run_infernal ) {

    /*
     * Rfam
     * doi: 10.1093/nar/gkx1038
     * url: https://rfam.xfam.org/
     *
     * Download Rfam to search with infernal.
     * Contains curated CMs/HMMs of ncRNA families.
     */
    process getRfam {

        label "download"
        label "small_task"
        time "4h"

        publishDir "${params.outdir}/downloads"

        output:
        set file("Rfam.cm"), file("Rfam.clanin") into rfam

        script:
        """
        wget -O Rfam.cm.gz "${params.rfam_url}"
        wget -O Rfam.clanin "${params.rfam_clanin_url}"
        gunzip Rfam.cm.gz
        """
    }
} else {

    rfam = Channel.empty()

}


if ( params.rfam_gomapping ) {
    Channel
        .value(file(params.rfam_gomapping, checkIfExists: true))
        .set { rfam2GO }
} else if ( run_infernal ) {

    process getRfam2GO {

        label "download"
        label "small_task"
        time "1h"

        publishDir "${params.outdir}/downloads"

        output:
        file "rfam2go" into rfam2GO

        script:
        """
        wget -O rfam2go "${params.rfam_gomapping_url}"
        """
    }
} else {
    rfam2GO = Channel.empty()
}


if ( params.pfam ) {
    Channel
        .fromPath( params.pfam, checkIfExists: true, type: "file")
        .collectFile(name: "pfam.stk", newLine: true)
        .set { pfamMSAs }
} else {

    /*
     * Download a number of Pfam domains and families associated with
     * TE activity.
     */
    process getPfam {

        label "download"
        label "medium_task"
        time "4h"

        publishDir "${params.outdir}/downloads"

        input:
        file "pfam_ids.txt" from pfamIds

        output:
        file "pfam.stk" into pfamMSAs

        script:
        """
        mkdir -p pfam
        cd pfam

        cat ../pfam_ids.txt \
        | xargs \
            -n 1 \
            -P "${task.cpus}" \
            -I "{}" \
            -- \
            wget \
              --no-check-certificate \
              -c \
              -O "{}.stk.gz" \
              "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{}/?annotation=alignment:full&download"

        gunzip *.gz
        cd ..
        cat pfam/*.stk > pfam.stk
        rm -rf -- pfam
        """
    }
}


if ( params.gypsydb ) {
    Channel
        .fromPath( params.gypsydb, checkIfExists: true, type: "file")
        .set { gypsydb }

} else {

    /*
     * GypsyDB
     *
     * doi: 10.1093/nar/gkq1061
     * url: http://www.gydb.org/index.php/Main_Page
     *
     * Curated msas and HMMs of many LTR families.
     * It apparently also contains some DNA elements.
     */
    process getGypsyDB {

        label "download"
        label "small_task"
        time "4h"

        publishDir "${params.outdir}/downloads"

        output:
        file "gypsydb/*" into gypsydb mode flatten

        script:
        """
        mkdir gypsydb
        wget --no-check-certificat -O gypsy.zip "${params.gypsydb_url}"
        unzip gypsy.zip
        mv GyDB*/alignments/*.sto ./gypsydb
        rm -rf -- GyDB*
        """
    }
}

process processGydb {

    label "posix"
    label "small_task"
    time "20m"

    input:
    file "stk/*" from gypsydb.collect()

    output:
    file "gypsydb.stk" into gypsyDBMSAs

    script:
    """
    touch gypsydb.stk
    for f in stk/*;
    do
      if grep -q "#=GF AC"
      then
        cat "\${f}" >> gypsydb.stk
      else
        # Strips directory and extension
        FILENAME=\$(basename \${f%.*})
        sed "1a #=GF AC \${FILENAME}" "\${f}" >> gypsydb.stk
      fi
    done
    """
}


Channel
    .fromPath(params.protein_families, checkIfExists: true, type: "file")
    .first()
    .set { proteinFamilies }


if ( params.mitefinder_profiles ) {
    Channel
        .fromPath(params.mitefinder_profiles, checkIfExists: true, type: "file")
        .first()
        .set { miteFinderProfiles }
} else {

    process getMiteFinderProfiles {

        label "mitefinder"
        label "small_task"
        time "20m"

        output:
        file "pattern_scoring.txt" into miteFinderProfiles

        script:
        """
        cp "\${MITEFINDER_PROFILE}" ./
        """
    }
}


genomes.into {
    genomes4RunTRNAScan;
    genomes4ChunkifyGenomes;
    genomes4RunStructRNAFinder;
    genomes4RunRNAmmer;
    genomes4RunOcculterCut;
    genomes4GetOcculterCutRegionFrequencies;
    genomes4GetOcculterCutGroupedRegionFrequencies;
    genomes4RunRepeatMaskerSpecies;
    genomes4RunRepeatMasker;
    genomes4RunRepeatModeler;
    genomes4GetMMSeqsGenomes;
    genomes4GetMMSeqsGenomeFastas;
    genomes4GetGtSuffixArrays;
    genomes4RunLtrHarvest;
    genomes4RunEAHelitron;
    genomes4RunMiteFinder;
    genomes4GetMiteFinderGFFs;
    genomes4GetMiteFinderFastas;
    genomes4GetSoftmaskedGenomes;
}


//
// Finding non-coding RNA
//


/*
 * tRNAScan-SE
 * doi: 10.1101/614032
 * url: http://lowelab.ucsc.edu/tRNAscan-SE/
 */
process runTRNAScan {

    label "trnascan"
    label "medium_task"
    time "5h"

    publishDir "${params.outdir}/${name}/noncoding"

    tag "${name}"

    input:
    set val(name), file(fasta) from genomes4RunTRNAScan

    output:
    set val(name),
        file("${name}_trnascan.txt"),
        file("${name}_trnascan_ss.txt") into tRNAScanResults

    set val(name), file("${name}_trnascan.fasta") into tRNAScanFasta

    file "${name}_trnascan_iso.txt"
    file "${name}_trnascan_stats.txt"
    file "${name}_trnascan.bed"

    script:
    """
    tRNAscan-SE \
      -E \
      -o "${name}_trnascan.txt" \
      -f "${name}_trnascan_ss.txt" \
      -s "${name}_trnascan_iso.txt" \
      -m "${name}_trnascan_stats.txt" \
      -b "${name}_trnascan.bed" \
      -a "${name}_trnascan.fasta" \
      --forceow \
      --log trna.log \
      --thread "${task.cpus}" \
      "${fasta}"
    """
}


/*
 * Convert trnascan-se output to nice GFF3.
 * Add predicted secondary structure and anticodons as attributes.
 */
process getTRNAScanGFF {

    label "gffpal"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    set val(name), file("ts.txt"), file("ss.txt") from tRNAScanResults

    output:
    set val(name),
        val("noncoding"),
        val("tRNAScan-SE"),
        val("trnascan"),
        file("trnascan.gff3") into tRNAScanGFF

    script:
    """
    gffpal trnascan2gff -o "trnascan.gff3" ts.txt ss.txt
    """
}


/*
 * Prepare the Rfam database to search with infernal.
 */
process pressRfam {

    label "infernal"
    label "small_task"
    time "2h"

    when:
    run_infernal

    input:
    set file("Rfam.cm"), file("Rfam.clanin") from rfam

    output:
    file "out" into pressedRfam

    script:
    """
    mkdir out
    cp -L Rfam.cm out/Rfam.cm
    cp -L Rfam.clanin out/Rfam.clanin
    cmpress -F out/Rfam.cm
    """
}


/*
 * Split genome fastas into multiple fastas for parallelisation.
 */
process chunkifyGenomes {

    label "python3"
    label "small_task"
    time "1h"

    tag "${name}"

    when:
    run_infernal

    input:
    set val(name),
        file("input.fasta") from genomes4ChunkifyGenomes

    output:
    set val(name),
        file("input.fasta"),
        file("out_*.fasta") into chunkifiedGenomes

    script:
    """
    chunk_genomes.py -n 32 --prefix "out_" input.fasta
    """
}


/*
 * Infernal
 * doi: 10.1093/bioinformatics/btt509
 * url: http://eddylab.org/infernal/
 *
 * Searches Rfam vs genome.
 * tblout2gff script comes from <https://github.com/nawrockie/jiffy-infernal-hmmer-scripts>
 */
process runInfernal {

    label "infernal"
    label "small_task"
    time "5h"

    tag "${name}"

    publishDir "${params.outdir}/${name}/noncoding", saveAs: exclude_unclean

    when:
    run_infernal

    input:
    set val(name),
        file("genome.fasta"),
        file("chunk.fasta") from chunkifiedGenomes
            .flatMap { n, g, cs -> cs.collect { c -> [n, g, c] } }

    file "rfam" from pressedRfam

    output:
    set val(name), file("cmscan.gff3.unclean") into infernalMatches
    file "${name}_rfam_cmscan.tblout"
    file "${name}_rfam_cmscan.out"

    script:
    """
    SIZE=\$(grep -v "^>" genome.fasta | sed 's/[[:space:]]//' | tr -d '\\n' | wc -c)
    TRSIZE=\$(perl -E "say \${SIZE} * 2 / 1e6")

    cmscan \
      --cpu 1 \
      -Z "\${TRSIZE}" \
      --cut_ga \
      --rfam \
      --nohmmonly \
      --tblout "${name}_rfam_cmscan.tblout" \
      --fmt 2 \
      --clanin rfam/Rfam.clanin \
      rfam/Rfam.cm \
      chunk.fasta \
    > "${name}_rfam_cmscan.out"

    grep -v ' = ' "${name}_rfam_cmscan.tblout" > filtered.txt

    infernal-tblout2gff.pl \
      --cmscan \
      --fmt2 \
      --all \
      -E "${params.infernal_max_evalue}" \
      filtered.txt \
    > "cmscan.gff3.unclean"

    rm filtered.txt
    """
}


process tidyInfernalMatches {

    label "gffpal"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    set val(name), file("infernal.gff3") from infernalMatches
    file "rfam2go" from rfam2GO

    output:
    set val(name), file("tidied.gff3") into tidiedInfernalMatches

    script:
    """
    awk -F '\\t' '
      BEGIN {OFS="\\t"}
      {\$9=gensub(":", "=", "g", \$9); print}
      ' infernal.gff3 \
    | tidy_infernal_gff.py \
      --go rfam2go \
      --best \
      --source "Rfam" \
      --type nucleotide_match \
      -o tidied.gff3 \
      -
    """
}


/*
 * Recombine infernal chunks into a single gff.
 */
process combineInfernal {

    label "genometools"
    label "small_task"
    time "1h"

    tag "${name}"

    when:
    run_infernal

    input:
    set val(name),
        file("unclean/*.gff3") from tidiedInfernalMatches
            .groupTuple(by: 0)

    output:
    set val(name),
        val("noncoding"),
        val("Rfam"),
        val("rfam_cmscan"),
        file("rfam_cmscan.gff3") into infernalGFF

    script:
    """
    mkdir clean

    for f in unclean/*.gff3
    do
      gt gff3 -tidy -sort "\${f}" > "clean/\$(basename \${f})"
    done

    gt merge -tidy clean/*.gff3 > "rfam_cmscan.gff3"
    """
}


/*
 * RNAmmer
 * doi: 10.1093/nar/gkm160
 * url: http://www.cbs.dtu.dk/services/RNAmmer/
 */
process runRNAmmer {

    label "rnammer"
    label "small_task"
    time "5h"

    publishDir "${params.outdir}/${name}/noncoding", saveAs: exclude_unclean

    tag "${name}"

    when:
    params.rnammer

    input:
    set val(name), file(fasta) from genomes4RunRNAmmer

    output:
    set val(name), file("rnammer.gff2.unclean") into rnammerResults
    file "${name}_rnammer.hmmreport"

    """
    rnammer \
      -S euk \
      -m lsu,ssu,tsu \
      -gff "rnammer.gff2.unclean" \
      -h "${name}_rnammer.hmmreport" \
      "${fasta}"
    """
}


/*
 * Convert gff2 format to gff3 and modify type labels.
 */
process getRNAmmerGFF {

    label "gffpal"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    set val(name), file("rnammer.gff2") from rnammerResults

    output:
    set val(name),
        val("noncoding"),
        val("RNAmmer"),
        val("rnammer"),
        file("rnammer.gff3") into rnammerGFF

    script:
    """
    gffpal rnammer2gff -k euk -o "rnammer.gff3" rnammer.gff2
    """
}


/*
 * OcculterCut
 * doi: 10.1093/gbe/evw121
 */
process runOcculterCut {

    label "occultercut"
    label "small_task"
    time "4h"
    tag "${name}"

    publishDir "${params.outdir}/${name}/noncoding"

    input:
    set val(name), file(fasta) from genomes4RunOcculterCut

    output:
    set val(name), file("${name}_occultercut_regions.gff3") into occulterCutRegions
    set val(name), file("${name}_occultercut_grouped_regions.gff3") optional true into occulterCutGroupedRegions
    file "${name}_occultercut.png" optional true
    file "${name}_occultercut_composition_gc.txt"
    file "${name}_occultercut_my_genome.txt"
    file "${name}_occultercut_nuc_frequencies.R*" optional true

    script:
    """
    OcculterCut -f "${fasta}"

    if [ -e plot.plt ]
    then
      sed -i '1i set terminal pngcairo size 900,600 enhanced font "Helvetica,20"' plot.plt
      sed -i '1a set output "plot.png"' plot.plt

      gnuplot plot.plt
      mv plot.png "${name}_occultercut.png"
    fi

    mv compositionGC.txt "${name}_occultercut_composition_gc.txt"
    mv regions.gff3 "${name}_occultercut_regions.gff3"
    mv myGenome.txt "${name}_occultercut_my_genome.txt"

    if [ -e groupedRegions.gff3 ]
    then
      mv groupedRegions.gff3 "${name}_occultercut_grouped_regions.gff3"
    fi

    find . -name "nuc_frequences.R*" -printf '%f\\0' \
    | xargs -I {} -0 -- mv '{}' "${name}_occultercut_{}"
    """
}


/*
 */
process getOcculterCutRegionFrequencies {

    label "gffpal"
    label "small_task"
    time "2h"

    tag "${name}"

    input:
    set val(name), file("in.gff"), file("in.fasta") from occulterCutRegions
        .combine(genomes4GetOcculterCutRegionFrequencies, by: 0)

    output:
    set val(name),
        val("occultercut_frequencies"),
        file("out.gff") into occulterCutRegionFrequencies

    script:
    """
    tidy_occultercut_regions.py \
      --source "OcculterCut" \
      --type "region" \
      -o out.gff \
      in.gff \
      in.fasta
    """
}


/*
 */
process getOcculterCutGroupedRegionFrequencies {

    label "gffpal"
    label "small_task"
    time "2h"

    tag "${name}"

    input:
    set val(name),
        file("in.gff"),
        file("in.fasta") from occulterCutGroupedRegions
            .combine(genomes4GetOcculterCutGroupedRegionFrequencies, by: 0)

    output:
    set val(name),
        val("occultercut_grouped_frequencies"),
        file("out.gff") into occulterCutGroupedRegionFrequencies

    script:
    """
    awk '
      BEGIN {OFS="\\t"; i=1}
      {print \$1, \$2, "region", \$4, \$5, \$6, \$7, \$8, "Name="\$3}
      ' in.gff \
    | tidy_occultercut_regions.py \
        --source "OcculterCut" \
        --type "region" \
        -o out.gff \
        - \
        in.fasta
    """
}


/*
 */
process tidyOcculterCutGFFs {

    label "genometools"
    label "small_task"
    time "1h"

    tag "${name} - ${suffix}"

    publishDir "${params.outdir}/${name}/noncoding"

    input:
    set val(name),
        val(suffix),
        file("in.gff") from occulterCutRegionFrequencies
            .mix(occulterCutGroupedRegionFrequencies)

    output:
    file "${name}_${suffix}.gff3"

    script:
    """
    gt gff3 \
      -tidy \
      -sort \
      in.gff \
    > "${name}_${suffix}.gff3"
    """
}


//
// Find repeats/TEs
//


/*
 * Prepare RepeatMasker database.
 */
process prepRepeatMaskerDB {

    label "repeatmasker"
    label "medium_task"
    time "6h"

    input:
    file "repbase.tar.gz" from repbase
    file "rm_meta.tar.gz" from rmMeta
    file "dfam38.0.h5.gz" from dfamH5
    file "Dfam.embl.gz" from dfamEMBLs
    file "RepeatPeps.lib" from rmRepeatPeps

    output:
    file "Libraries" into rmlib

    script:
    """
    cp -r "\${RMASK_PREFIX}/Libraries" Libraries

    # These tar archives unpack to a folder "Libraries"
    if [ -e "repbase.tar.gz" ]
    then
        tar zxf "repbase.tar.gz"
    fi

    if [ -e "rm_meta.tar.gz" ]
    then
        tar zxf "rm_meta.tar.gz"
    fi

    # Repeat peps gets distributed with the repeat masker executables.
    # Not available elsewhere.
    cp -L RepeatPeps.lib Libraries
    cp -L Dfam.embl.gz Libraries
    mkdir -p Libraries/famdb
    cp -L dfam38.0.h5.gz Libraries/famdb

    cd Libraries

    gunzip --force famdb/dfam38.0.h5.gz
    # DFAM consensus filenames hard-coded into the script, so we move it.
    mv Dfam.embl.gz DfamConsensus.embl.gz
    gunzip --force DfamConsensus.embl.gz

    # I think this comes with meta.
    if [ -e "taxonomy.dat.gz" ]
    then
      gunzip --force taxonomy.dat.gz
    fi

    # This basically concatenates the different blastable (i.e. not hmm) databases together.
    perl \
      -I "\${RMASK_PREFIX}" \
      -e "use LibraryUtils; LibraryUtils::rebuildMainLibrary( \\"../Libraries\\" );"
    
    sed -i '/^\$/d' RepeatMaskerLib.embl
    buildRMLibFromEMBL.pl RepeatMaskerLib.embl > RepeatMasker.lib

    makeblastdb -dbtype nucl -in RepeatMasker.lib
    makeblastdb -dbtype prot -in RepeatPeps.lib

    # As far as I can tell, nothing needs to be done with dfam hmms?
    """
}


/*
 */
process filterScaffoldLength {

    label "posix"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    set val(name), file("in.fasta") from genomes4RunRepeatModeler

    output:
    set val(name), file("filtered.fasta") into filteredGenomes4RunRepeatModeler

    script:
    """
    fasta_to_tsv.sh < in.fasta \
    | awk -v size="${params.repeatmodeler_min_len}" 'length(\$2) >= size' \
    | tsv_to_fasta.sh \
    > filtered.fasta
    """
}


/*
 * RepeatModeler
 * url: http://www.repeatmasker.org/RepeatModeler/
 *
 * Note that if it doesn't find any families, then
 * RepeatModeler won't output any files. It has to be optional.
 */
process runRepeatModeler {

    label "repeatmasker"
    label "big_task"
    time "1d"

    errorStrategy 'retry'
    maxRetries 3

    tag "${name}"
    publishDir "${params.outdir}/${name}/tes"

    input:
    set val(name), file(fasta) from filteredGenomes4RunRepeatModeler
    file "rmlib" from rmlib

    output:
    set val(name),
        file("${name}_repeatmodeler.stk") optional true into repeatModelerStks
    file "${name}_repeatmodeler_consensus.fasta" optional true

    script:
    // Task uses n+1 threads.
    def ncpu = task.cpus == 1 ? 1 : task.cpus - 1
    """
    # repeatmasker modifies the content of rmlib.
    # I'm not sure if repeatmodeler will too, but just to be safe.
    # We copy it to avoid messing up checkpointing.

    # Where will this be placed?
    BuildDatabase -name "${name}" -engine ncbi "${fasta}"

    # Can we split this over scaffolds?
    RepeatModeler -engine ncbi -threads "${ncpu}" -database "${name}" >& run.out

    if [ -e "${name}-families.fa" ]
    then
      mv "${name}-families.fa" "${name}_repeatmodeler_consensus.fasta"
    fi

    if [ -e "${name}-families.stk" ]
    then
      mv "${name}-families.stk" "${name}_repeatmodeler.stk"
    fi
    """
}


repeatModelerStks.into {
    repeatModelerStks4GetFasta;
    repeatModelerStks4GetGFF;
}


/*
 * Get the ungapped fasta-formatted sequences from repeatmodeler so
 * we can cluster them.
 */
process getRepeatModelerFasta {

    label "python3"
    label "small_task"
    time "1h"

    tag "${name}"
    publishDir "${params.outdir}/${name}/tes"

    input:
    set val(name), file("in.stk") from repeatModelerStks4GetFasta

    output:
    set val(name), file("${name}_repeatmodeler.fasta") into repeatModelerFasta

    script:
    """
    stk2fasta.py "in.stk" \
    | sed "s/^>/>${name}:/" \
    > "${name}_repeatmodeler.fasta"
    """
}


/*
 * Convert the RepeatModeler STK to a GFF3 (as best as we can).
 */
process getRepeatModelerGFF {

    label "gffpal"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    set val(name), file("in.stk") from repeatModelerStks4GetGFF

    output:
    set val(name),
        val("tes"),
        val("RepeatModeler"),
        val("repeatmodeler"),
        file("repeatmodeler.gff3") into repeatModelerGFF

    script:
    """
    repeatmodeler2gff.py \
      -s RepeatModeler \
      -t repeat_region \
      -o repeatmodeler.gff3 \
      in.stk
    """
}


/*
 * Index genomes for MMSeqs
 */
process getMMSeqsGenomes {

    label "mmseqs"
    label "small_task"
    time "2h"

    tag "${name}"

    input:
    set val(name), file(fasta) from genomes4GetMMSeqsGenomes

    output:
    set val(name), file("genome") into mmseqsGenomes

    script:
    """
    mkdir genome orfs translated_orfs tmp

    mmseqs createdb "${fasta}" genome/db #--dont-split-seq-by-len

    # This command seems to be the secret sauce to getting profile-vs-genome
    # searches to work.
    mmseqs createindex genome/db tmp --threads "${task.cpus}" --search-type 2
    rm -rf -- tmp
    """
}


pfamMSAs.map { ["Pfam", it] }
    .mix(
        gypsyDBMSAs.map { ["GyDb", it] },
        proteinFamilies.map { ["pante_protein_families", it] },
    )
    .into { msas4GetAttributes; msas4GetProfiles }


/*
 * Extract information from the stockholm alignments to use for GFF3
 * attributes later.
 */
process getMSAAttributes {

    label "gffpal"
    label "small_task"
    time "1h"

    tag "${db}"

    input:
    set val(db), file("msas.stk") from msas4GetAttributes

    output:
    set val(db), file("msas_attributes.tsv") into msaAttributes

    script:
    """
    pfam2tab.py msas.stk > msas_attributes.tsv
    """
}


/*
 * Get MSA PSSM profile for each family.
 */
process getMMSeqsProfiles {

    label "mmseqs"
    label "medium_task"
    time "4h"

    tag "${db}"

    input:
    set val(db), file("msas.stk") from msas4GetProfiles

    output:
    set val(db), file("profiles") into mmseqsProfiles

    script:
    """
    mkdir msas profiles tmp
    mmseqs convertmsa msas.stk msas/db --identifier-field 1

    mmseqs msa2profile \
      msas/db \
      profiles/db \
      --match-mode 1 \
      --match-ratio 0.5 \
      --threads "${task.cpus}"

    #mmseqs createindex profiles/db tmp -k 6 -s 7 --threads "${task.cpus}"

    rm -rf -- tmp msas
    """
}


/*
 * MMSeqs2
 * doi: https://www.nature.com/articles/nbt.3988
 * url: https://github.com/soedinglab/MMseqs2
 *
 * Search the profiles against the genome orfs.
 * We could potentially look at "enriching" the profiles using extracted orfs.
 * I've tried it, it works, but it requires re-aligning the matching sequences
 * to the profiles outside of mmseqs, which is cumbersome.
 */
process searchProfilesVsGenomes {

    label "mmseqs"
    label "medium_task"
    time "4h"

    tag "${name} - ${db}"
    publishDir "${params.outdir}/${name}/tes"

    input:
    set val(name),
        file("genome"),
        val(db),
        file("profiles") from mmseqsGenomes.combine(mmseqsProfiles)

    output:
    set val(db),
        val(name),
        file("${name}_${file_db}_search.tsv") into mmseqsGenomeMatches

    script:
    file_db = db.toLowerCase()

    """
    mkdir search tmp

    mmseqs search \
      profiles/db \
      genome/db \
      search/db \
      tmp \
      --threads "${task.cpus}" \
      -e "${params.mmseqs_max_evalue}" \
      -s 7.5 \
      --search-type 2 \
      --num-iterations 1 \
      --min-length 10 \
      --mask 0 \
      --orf-start-mode 1 \
      --realign \
      -a

    mmseqs convertalis \
      profiles/db \
      genome/db \
      search/db \
      search_tmp.tsv \
      --threads "${task.cpus}" \
      --format-mode 0 \
      --format-output 'target,query,tstart,tend,tlen,qstart,qend,qlen,evalue,gapopen,pident,alnlen,raw,bits,cigar,mismatch,qcov,tcov'

    sort -k1,1 -k3,3n -k4,4n -k2,2 search_tmp.tsv > "${name}_${file_db}_search.tsv"
    sed -i '1i #target\tquery\ttstart\ttend\ttlen\tqstart\tqend\tqlen\tevalue\tgapopen\tpident\talnlen\traw\tbits\tcigar\tmismatch\tqcov\ttcov' "${name}_${file_db}_search.tsv"
    """
}


/*
 * Convert MMSeqs matches to a gff3 file.
 */
process getMMSeqsGenomeGFFs {

    label "gffpal"
    label "small_task"
    time "1h"

    tag "${name} - ${db}"

    input:
    set val(db),
        val(name),
        file("search.tsv"),
        file("attributes.tsv") from mmseqsGenomeMatches
            .combine(msaAttributes, by: 0)

    output:
    set val(name),
        val(db),
        file("search.gff3") into mmseqsGenomeGFFUnclean

    script:
    """
    mmseqs2gff.py \
      -s "${db}" \
      -a attributes.tsv \
      -o "search.gff3" \
      search.tsv
    """
}


/*
 * Get the sequences and make the GFF3 compliant.
 */
process getMMSeqsGenomeFastas {

    label "genometools"
    label "small_task"
    time "1h"

    tag "${name} - ${db}"

    publishDir "${params.outdir}/${name}/tes", saveAs: exclude_unclean

    input:
    set val(name),
        val(db),
        file("input.gff3"),
        file("genome.fasta") from mmseqsGenomeGFFUnclean
            .combine(genomes4GetMMSeqsGenomeFastas, by: 0)

    output:
    set val(name),
        val("tes"),
        val(db),
        val("${file_db}_search"),
        file("search.gff3.unclean") into mmseqsGenomeGFF

    set val(name),
        file("${name}_${file_db}_search.fasta") into mmseqsGenomeFasta

    script:
    file_db = db.toLowerCase()

    """
    gt gff3 \
      -tidy \
      -sort \
      input.gff3 \
    > search.gff3.unclean

    gt extractfeat \
      -type nucleotide_to_protein_match \
      -matchdescstart \
      -seqid \
      -coords \
      -seqfile genome.fasta \
      search.gff3.unclean \
    | fix_fasta_names.sh "${name}" \
    > "${name}_${file_db}_search.fasta"
    """
}


/*
 * Both ltrharvest and ltrdigest use these suffix arrays
 */
process getGtSuffixArrays {

    label "genometools"
    label "small_task"
    time "2h"

    tag "${name}"

    input:
    set val(name), file(fasta) from genomes4GetGtSuffixArrays

    output:
    set val(name),
        file(fasta),
        file("${fasta}*") into gtSuffixArrays

    script:
    """
    # Create the suffix arrays
    gt suffixerator \
      -db "${fasta}" \
      -indexname "${fasta}" \
      -lossless \
      -tis \
      -suf \
      -lcp \
      -des \
      -ssp \
      -sds \
      -dna
    """
}

gtSuffixArrays.into {
    gtSuffixArrays4RunLtrHarvest;
    gtSuffixArrays4RunLtrDigest;
}


/*
 * LTRHarvest
 * doi: 10.1186/1471-2105-9-18
 *
 * We may need to split this into two commands
 * 1 to find recent ones, and 1 to find old ones.
 * The options '-motif tgca' is usually used for canonical ones.
 */
process runLtrHarvest {

    label "genometools"
    label "small_task"
    time "4h"

    tag "${name}"

    input:
    set val(name), file(fasta), file("*") from gtSuffixArrays4RunLtrHarvest

    output:
    set val(name), file("ltrharvest.gff3") into ltrHarvestResults

    """
    gt ltrharvest \
      -seqids yes \
      -similar "${params.ltrharvest_similar}" \
      -vic "${params.ltrharvest_vic}" \
      -seed "${params.ltrharvest_seed}" \
      -minlenltr "${params.ltrharvest_minlenltr}" \
      -maxlenltr "${params.ltrharvest_maxlenltr}" \
      -mintsd "${params.ltrharvest_mintsd}" \
      -maxtsd "${params.ltrharvest_maxtsd}" \
      -index "${fasta}" \
      -gff3 "ltrharvest_tmp.gff3"

    gt gff3 \
      -tidy \
      -sort \
      -retainids \
      "ltrharvest_tmp.gff3" \
    > "ltrharvest.gff3"
    """
}


/*
 * LTRDigest
 * doi: 10.1093/nar/gkp759
 *
 * TODO: filter out incomplete LTRs?
 * This should remove many false positives, but might might exclude things.
 *
 * It seems like the HMM options don't work anymore.
 */
process runLtrDigest {

    label "genometools"
    label "small_task"
    time "4h"

    tag "${name}"
    publishDir "${params.outdir}/${name}/tes", saveAs: exclude_unclean

    input:
    set val(name),
        file(fasta),
        file("*"),
        file("${name}_ltrharvest.gff3"),
        file("trna.fasta") from gtSuffixArrays4RunLtrDigest
            .combine(ltrHarvestResults, by: 0)
            .combine(tRNAScanFasta, by: 0)

    output:
    set val(name), file("${name}_ltrdigest.fasta") into ltrDigestFasta
    set val(name),
        val("tes"),
        val("KEEP"),
        val("ltrdigest"),
        file("ltrdigest.gff3.unclean") into ltrDigestGFF

    file "${name}_ltrdigest_complete.fas"
    file "${name}_ltrdigest_3ltr.fas"
    file "${name}_ltrdigest_5ltr.fas"
    file "${name}_ltrdigest_conditions.csv"
    file "${name}_ltrdigest_tabout.csv"

    script:
    """
    mkdir -p tmp
    TMPDIR="\${PWD}/tmp" \
    gt -j "${task.cpus}" ltrdigest \
      -trnas trna.fasta \
      -outfileprefix "${name}_ltrdigest" \
      -matchdescstart \
      -seqfile "${fasta}" \
      "${name}_ltrharvest.gff3" \
    > "ltrdigest_tmp.gff3"

    gt gff3 \
      -tidy \
      -sort \
      -retainids \
      "ltrdigest_tmp.gff3" \
    > "ltrdigest.gff3.unclean"

    gt extractfeat \
      -type repeat_region \
      -matchdescstart \
      -seqid \
      -coords \
      -seqfile "${fasta}" \
      "ltrdigest.gff3.unclean" \
    | fix_fasta_names.sh "${name}" \
    > "${name}_ltrdigest.fasta"

    rm -rf -- tmp
    """
}


/*
 * EAHelitron
 * doi: 10.1186/s12859-019-2945-8
 * url: https://github.com/dontkme/EAHelitron
 *
 * TODO: compare performance with HelitronScanner.
 * Fungal helitrons might have non-standard motifs.
 * Parallel version doesn't support all arguments that the single core version does.
 * See doi: 10.1186/s13100-016-0083-7
 */
process runEAHelitron {

    label "eahelitron"
    label "small_task"
    time "6h"

    tag "${name}"
    publishDir "${params.outdir}/${name}/tes", saveAs: exclude_unclean

    input:
    set val(name), file(fasta) from genomes4RunEAHelitron

    output:
    set val(name), file("${name}_eahelitron.fasta") into eaHelitronFasta
    set val(name), file("eahelitron.gff3.unclean") into eaHelitronUncleanGFF
    file "${name}_eahelitron.5.fa"
    file "${name}_eahelitron.3.txt"
    file "${name}_eahelitron.5.txt"
    file "${name}_eahelitron.u*.fa"
    file "${name}_eahelitron.d*.fa"
    file "${name}_eahelitron.bed"
    file "${name}_eahelitron.len.txt"

    script:
    """
    EAHelitron \
      -o "${name}_eahelitron" \
      -r "${params.eahelitron_three_prime_fuzzy_level}" \
      -u "${params.eahelitron_upstream_length}" \
      -d "${params.eahelitron_downstream_length}" \
      "${fasta}"

    awk -v name="${name}" '
      /^>/ {
        pos=gensub(/^>\\S+\\s+([^:]+):([0-9]+)\\.\\.([0-9]+)\$/, "\\\\1:\\\\2-\\\\3", "g", \$0);
        \$0 = ">" name ":" pos;
      }
      { print }
    ' < "${name}_eahelitron.5.fa" \
    > "${name}_eahelitron.fasta"

    mv "${name}_eahelitron.gff3" "eahelitron.gff3.unclean"
    """
}


process filterEAHelitronGFF {

    label "gffpal"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    set val(name), file("in.gff") from eaHelitronUncleanGFF

    output:
    set val(name),
        val("tes"),
        val("EAhelitron"),
        val("eahelitron"),
        file("eahelitron.gff3") into eaHelitronGFF

    script:
    """
    filter_eahelitron_gff.py -o "eahelitron.gff3" in.gff
    """
}


/*
 * MITEfinderII
 * doi: 10.1186/s12920-018-0418-y
 * url: https://github.com/screamer/miteFinder
 *
 * TODO: Compare performance with MiteTracker if can resolve the
 * lousy install process.
 */
process runMiteFinder {
    label "mitefinder"
    label "small_task"
    time "6h"

    tag "${name}"

    input:
    set val(name), file(genome) from genomes4RunMiteFinder
    file "pattern_scoring.txt" from miteFinderProfiles

    output:
    set val(name), file("mitefinder.fasta") into miteFinderResults

    script:
    """
    miteFinder \
      -input "${genome}" \
      -output "mitefinder.fasta" \
      -pattern_scoring pattern_scoring.txt \
      -threshold "${params.mitefinder_threshold}"
    """
}


/*
 * Construct a gff from the mitefinder fasta descriptions.
 */
process getMiteFInderGFFs {

    label "gffpal"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    set val(name),
        file("in.fasta"),
        file("genome.fasta") from miteFinderResults
            .combine(genomes4GetMiteFinderGFFs, by: 0)

    output:
    set val(name),
        file("out.gff3") into miteFinderGFF

    script:
    """
    mitefinder2gff.py genome.fasta in.fasta > out.gff3
    """
}


miteFinderGFF
    .tap { miteFinderGFF4GetMiteFinderFasta }
    .map { n, g -> [n, "tes", "MiteFinderII", "mitefinder", g] }
    .set { miteFinderGFF4TidyGFF }


/*
 * Make the mitefinder gffs compliant, and extract new sequences
 * with new names.
 */
process getMiteFinderFastas {

    label "genometools"
    label "small_task"
    time "1h"

    tag "${name}"
    publishDir "${params.outdir}/${name}/tes"

    input:
    set val(name),
        file("in.gff3"),
        file("genome.fasta") from miteFinderGFF4GetMiteFinderFasta
            .combine(genomes4GetMiteFinderFastas, by: 0)

    output:
    set val(name),
        file("${name}_mitefinder.fasta") into miteFinderFasta

    script:
    """
    gt gff3 \
      -tidy \
      -sort \
      "in.gff3" \
    | gt extractfeat \
      -type repeat_region \
      -matchdescstart \
      -seqid \
      -coords \
      -seqfile "genome.fasta" \
      - \
    | fix_fasta_names.sh "${name}" \
    > "${name}_mitefinder.fasta"
    """
}


/*
 * This just concatenates all fasta files, removes sequences with
 * exactly the same name, which should all have the same sequence
 * because they are named by the genomic coords that they were taken from.
 *
 * I'm leaving out LTRdigest matches because they are usually very big,
 * So would tend to merge tandem fragments.
 * TODO: filter out overlapping matches to have a single one per locus.
 */
process combineTEFastas {

    label "python3"
    label "small_task"
    time "2h"

    publishDir "${params.outdir}/pantes"

    input:
    file "in/*.fasta" from repeatModelerFasta
        .mix(
            eaHelitronFasta,
            miteFinderFasta,
            mmseqsGenomeFasta,
        )
        .map { n, f -> f }
        .collect()

    output:
    file "combined_tes.fasta" into combinedTEFasta

    script:
    """
    cat in/*.fasta \
    | exclude_contained_subseqs.py \
        --coverage 0.95 \
        -o combined_tes.fasta \
        -
    """
}


/*
 * vsearch
 * doi: 10.7717/peerj.2584
 * url: https://github.com/torognes/vsearch
 *
 * Cluster to find naive families.
 * After this point, we could filter by frequency per-genome and within the
 * population. This requires that we de-duplicate loci first.
 */
process clusterTEFastas {

    label "vsearch"
    label "medium_task"
    time "1d"

    publishDir "${params.outdir}/pantes"

    input:
    file "combined.fasta" from combinedTEFasta

    output:
    set file("clusters.tsv"), file("clusters") into clusteredTEFasta

    script:
    """
    mkdir -p clusters
    vsearch \
      --threads "${task.cpus}" \
      --cluster_fast "combined.fasta" \
      --id 0.90 \
      --weak_id 0.7 \
      --iddef 0 \
      --qmask dust \
      --uc clusters.tsv \
      --strand "both" \
      --clusters "clusters/fam"
    """
}


/*
 * Filter clusters based on frequencies.
 */
process filterTEClusters {

    label "python3"
    label "small_task"
    time "4h"

    publishDir "${params.outdir}/pantes"

    input:
    set file("clusters.tsv"), file("clusters") from clusteredTEFasta

    output:
    file "filtered_clusters_included" into filteredTEClusters
    file "filtered_clusters_excluded"
    file "filtered_clusters.tsv"

    script:
    """
    filter_clusters.py \
      --outsel filtered_clusters_included \
      --outexcl filtered_clusters_excluded \
      --outtab filtered_clusters.tsv \
      --min-intra "${params.min_intra_frequency}" \
      --min-inter "${params.min_inter_proportion}" \
      clusters.tsv \
      clusters/*
    """
}


/*
 * DECIPHER
 * doi: 10.1186/s12859-015-0749-z
 * url: http://www2.decipher.codes/AlignSequences.html
 *
 * Multiple sequence alignment of the clusters.
 *
 * To expand clusters to include more distant matches, we could
 * convert these MSAs to HMMs and search HMM vs consensus sequences
 * using nhmmer?
 */
process getClusterMSAs {

    label "decipher"
    label "big_task"
    time "1d"

    publishDir "${params.outdir}/pantes"

    input:
    file "clusters" from filteredTEClusters

    output:
    file "msas" into clusterMSA

    script:
    """
    mkdir -p msas

    find clusters/ -name "fam*" -printf '%f\\0' \
    | xargs -0 -P "${task.cpus}" -I {} -- \
        run_decipher.R \
          --infile "clusters/{}" \
          --outfile "msas/{}"
    """
}


clusterMSA.into {
    clusterMSA4GetClusterMSAConsensus;
    clusterMSA4GetClusterMSAStockholm;
}


/*
 * Find consensus sequences of the MSAs.
 *
 * This just takes a majority-rules approach to consensus finding.
 * Gaps are ignored in the calculation so no gaps should appear in the consensus.
 * Ties are broken by selecting the first match in the following order
 * N, A, C, G, T
 */
process getClusterMSAConsensus {

    label "decipher"
    label "medium_task"
    time "6h"

    publishDir "${params.outdir}/pantes"

    input:
    file "msas" from clusterMSA4GetClusterMSAConsensus

    output:
    file "families_consensi.fasta" into clusterMSAConsensus

    script:
    """
    mkdir consensi

    find msas/ -name "fam*" -printf '%f\\0' \
    | xargs -0 -P "${task.cpus}" -I {} -- \
        get_consensus.R \
          --infile "msas/{}" \
          --outfile "consensi/{}"

    cat consensi/* > families_consensi.fasta
    rm -rf -- consensi
    """
}


/*
 * Convert the fasta multiple sequence alignments to a single stockholm MSA
 * file.
 */
process getClusterMSAStockholm {

    label "python3"
    label "small_task"
    time "2h"

    publishDir "${params.outdir}/pantes"

    input:
    file "msas" from clusterMSA4GetClusterMSAStockholm

    output:
    file "families.stk" into clusterMSAStockholm

    script:
    """
    fasta_aln2stk.py -o families.stk msas/*
    """
}


/*
 * RepeatClassifier
 * url: http://www.repeatmasker.org/RepeatModeler/
 *
 * Classifies the custom repeat library by matches to the consensus
 * sequences. Also adds annotations to the stockholm, but AFAIK this
 * is just transferred from the consensi.
 */
process runRepeatClassifier {

    label "repeatmasker"
    label "big_task"
    time "12h"

    publishDir "${params.outdir}/pantes"

    input:
    file "families_consensi.fasta" from clusterMSAConsensus
    file "families.stk" from clusterMSAStockholm
    file "rmlib" from rmlib

    output:
    file "families_classified_consensi.fasta" into repeatClassifierResults
    file "families_classified.stk"

    script:
    """
    RepeatClassifier \
      -consensi families_consensi.fasta \
      -stockholm families.stk

    mv families_consensi.fasta.classified families_classified_consensi.fasta
    mv families-classified.stk families_classified.stk
    """
}


/*
 * RepeatMasker
 * url: http://www.repeatmasker.org/RMDownload.html
 *
 * This uses pre-existing species information for repeatmasker.
 * Only run if user provides a species to use.
 */
process runRepeatMaskerSpecies {

    label "repeatmasker"
    label "medium_task"
    time "1d"

    tag "${name}"
    publishDir "${params.outdir}/${name}/tes"

    when:
    params.rm_species

    input:
    set val(name), file(fasta) from genomes4RunRepeatMaskerSpecies
    file "rmlib" from rmlib

    output:
    set val(name), file("${name}_repeatmasker_species.txt") into repeatMaskerSpeciesResults
    file "${name}_repeatmasker_species_align.txt"
    file "${name}_repeatmasker_species_masked.fasta"
    file "${name}_repeatmasker_species_repeat_densities.txt"

    script:
    """
    # repeatmasker modifies the content of rmlib.
    # We copy it to avoid messing up checkpointing.
    cp -rL rmlib rmlib_tmp

    RepeatMasker \
      -e ncbi \
      -pa "${task.cpus}" \
      -species "${params.rm_species}" \
      -libdir "./rmlib_tmp" \
      -xsmall \
      -gff \
      -alignments \
      -excln \
      "${fasta}"

    mv "${fasta.name}.align" "${name}_repeatmasker_species_align.txt"
    mv "${fasta.name}.masked" "${name}_repeatmasker_species_masked.fasta"
    mv "${fasta.name}.tbl" "${name}_repeatmasker_species_repeat_densities.txt"
    mv "${fasta.name}.out" "${name}_repeatmasker_species.txt"
    rm -rf -- ./rmlib_tmp
    """
}


/*
 * RepeatMasker
 * url: http://www.repeatmasker.org/RMDownload.html
 */
process runRepeatMasker {

    label "repeatmasker"
    label "medium_task"
    time "1d"

    tag "${name}"
    publishDir "${params.outdir}/${name}/tes"

    input:
    set val(name), file(fasta) from genomes4RunRepeatMasker
    file "families_consensi.fasta" from repeatClassifierResults
    file "rmlib" from rmlib

    output:
    set val(name), file("${name}_repeatmasker.txt") into repeatMaskerResults
    file "${name}_repeatmasker_align.txt"
    file "${name}_repeatmasker_masked.fasta"
    file "${name}_repeatmasker_repeat_densities.txt"

    script:
    """
    # repeatmasker modifies the content of rmlib.
    # We copy it to avoid messing up checkpointing.
    cp -rL rmlib rmlib_tmp

    RepeatMasker \
      -e ncbi \
      -pa "${task.cpus}" \
      -lib families_consensi.fasta \
      -libdir "./rmlib_tmp" \
      -xsmall \
      -alignments \
      -excln \
      "${fasta}"

    mv "${fasta.name}.align" "${name}_repeatmasker_align.txt"
    mv "${fasta.name}.masked" "${name}_repeatmasker_masked.fasta"
    mv "${fasta.name}.tbl" "${name}_repeatmasker_repeat_densities.txt"
    mv "${fasta.name}.out" "${name}_repeatmasker.txt"
    rm -rf -- ./rmlib_tmp
    """
}


process getRepeatMaskerGFF {

    label "gffpal"
    label "small_task"
    time "2h"

    tag "${name} - ${analysis}"

    input:
    set val(name),
        val(analysis),
        file("rm.out") from repeatMaskerResults
            .map { n, o -> [n, "repeatmasker", o] }
            .mix( repeatMaskerSpeciesResults
                    .map {n, o -> [n, "repeatmasker_species", o ]} )

    output:
    set val(name),
        val("tes"),
        val("RepeatMasker"),
        val(analysis),
        file("rm.gff3") into repeatMaskerGFF

    script:
    """
    rmout2gff3.py -o "rm.gff3" rm.out
    """
}


process tidyGFFs {

    label "genometools"
    label "small_task"
    time "1h"

    publishDir "${params.outdir}/${name}/${folder}"

    tag "${name} - ${analysis}"

    input:
    set val(name),
        val(folder),
        val(source),
        val(analysis),
        file("in.gff") from tRNAScanGFF
            .mix(
                infernalGFF,
                rnammerGFF,
                repeatModelerGFF,
                mmseqsGenomeGFF,
                ltrDigestGFF,
                eaHelitronGFF,
                miteFinderGFF4TidyGFF,
                repeatMaskerGFF,
            )

    output:
    set val(name),
        val(folder),
        file("${name}_${analysis}.gff3") into tidiedGFFs

    script:
    setsource = source == "KEEP" ? "" : "-setsource ${source}"

    """
    grep -v "^#" in.gff \
    | gt gff3 \
        -tidy \
        -sort \
        ${setsource} \
        -retainids \
        - \
    > "${name}_${analysis}.gff3"
    """
}


/*
 * This just combines all results into some final GFFs
 */
process combineGFFs {

    label "genometools"
    label "small_task"
    time "2h"

    publishDir "${params.outdir}/${name}/final"
    tag "${name}"

    input:
    set val(name),
        val(suffix),
        file("gffs/*.gff3") from tidiedGFFs
            .flatMap { n, f, g -> [[n, "", g], [n, "_${f}", g]] }
            .groupTuple(by: [0, 1])

    output:
    set val(name),
        val(suffix),
        file("${name}_pante${suffix}.gff3") into combinedGFF

    script:
    """
    gt merge -tidy gffs/*.gff3 > "${name}_pante${suffix}.gff3"
    """
}


/*
 * Create a soft-masked fasta file given our repeat annotations.
 */
process getSoftmaskedGenomes {

    label "bedtools"
    label "small_task"
    time "1h"

    publishDir "${params.outdir}/${name}/final"
    tag "${name}"

    input:
    set val(name),
        file("genome.fasta"),
        file("repeats.gff3") from genomes4GetSoftmaskedGenomes
            .combine(
                combinedGFF
                    .filter { n, s, f -> s == "" }
                    .map { n, s, f -> [n, f] },
                by: 0
            )

    output:
    file "${name}_softmasked.fasta"

    script:
    """
    fasta_to_tsv.sh \
    < genome.fasta \
    | awk 'BEGIN {OFS="\\t"} {print \$1, 0, length(\$2)}' \
    > genome.bed

    bedtools intersect \
      -a repeats.gff3 \
      -b genome.bed \
    > bounded.gff3


    bedtools maskfasta \
      -fi genome.fasta \
      -bed bounded.gff3 \
      -fo "${name}_softmasked.fasta" \
      -soft

    rm -f bounded.gff3 genome.bed
    """
}
