#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # pante


    Example:
    nextflow run -resume -profile dev,standard main.nf \
      --genomes "A_fumigatus_A1163_chr1.fasta" \
      --repbase "containers/rmlib" \
      --rmspecies "fungi" \
      --rnammer


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
params.repbase = false
params.rmspecies = "fungi"
params.helitronscanner_heads = "$baseDir/data/helitronscanner_head.lcvs"
params.helitronscanner_tails = "$baseDir/data/helitronscanner_tail.lcvs"

params.mitefinder_profiles = false

params.structrnafinder = false
params.rfam = false
params.rfam_url = "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
params.rnammer = false

params.ltr_hmms = false
params.gypsydb_url = "http://gydb.org/gydbModules/collection/collection/db/GyDB_collection.zip"
params.pfam_ltr_ids = "$baseDir/data/pfam_ltr_domains.txt"


if ( params.genomes ) {
    genomes = Channel
        .fromPath(params.genomes, checkIfExists: true, type: "file")
        .map { file -> [file.baseName, file] }
} else {
    log.info "Hey I need some genomes to assess please."
    exit 1
}


if ( params.repbase ) {
    repbase = Channel.fromPath(params.repbase, checkIfExists: true).first()
} else {
    log.info "Sorry for now we need repbase"
    exit 1
}


if ( params.pfam_ltr_ids ) {
    Channel
        .fromPath(params.pfam_ltr_ids, checkIfExists: true, type: "file")
        .first()
        .set { pfamLtrIds }
} else {
    pfamLtrIds = Channel.empty()
}


if ( params.rfam ) {
    Channel
        .fromPath( params.rfam, checkIfExists: true, type: "file")
        .first()
        .set { rfam }
} else if ( params.structrnafinder ) {

    process getRfam {

        label "download"
        label "small_task"

        publishDir "${params.outdir}/downloads"

        output:
        file "Rfam.cm" into rfam

        script:
        """
        wget -O Rfam.cm.gz "${params.rfam_url}"
        gunzip Rfam.cm.gz
        """
    }
} else {
    rfam = Channel.empty()
}


if ( params.ltr_hmms ) {
    Channel
        .fromPath( params.ltr_hmms, checkIfExists: true, type: "file")
        .collect()
        .set { ltrHMMs }
} else {

    process getGypsyDB {

        label "download"
        label "small_task"

        publishDir "${params.outdir}/downloads"

        output:
        file "gypsydb/*" into gypsydbHMMs

        script:
        """
        wget -O gypsy.zip "${params.gypsydb_url}"
        unzip gypsy.zip
        mv GyDB*/profiles ./gypsydb
        rm -rf -- GyDB*
        """
    }


    process getPfamHmms {

        label "download"
        label "small_task"

        publishDir "${params.outdir}/downloads"

        input:
        file "pfam_ids.txt" from pfamLtrIds

        output:
        file "pfam_ltrs/*" into pfamLtrHMMs mode flatten

        script:
        """
        mkdir -p pfam_ltrs

          cat pfam_ids.txt \
        | xargs \
            -n 1 \
            -P "${task.cpus}" \
            -I "{}" \
            -- \
            wget --no-check-certificate -P pfam_ltrs "https://pfam.xfam.org/family/{}/hmm"
        """
    }

    gypsydbHMMs.mix(pfamLtrHMMs).collect().set { ltrHMMs }
}


if ( params.mitefinder_profiles ) {
    Channel
        .fromPath(params.mitefinder_profiles, checkIfExists: true, type: "file")
        .first()
        .set { miteFinderProfiles }
} else {

    process getMiteFinderProfiles {

        label "mitefinder"
        label "small_task"

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
    genomes4RunStructRNAFinder;
    genomes4RunRNAmmer;
    genomes4RunOcculterCut;
    genomes4RunRepeatMaskerRepbase;
    genomes4RunRepeatModeller;
    genomes4RunLtrHarvest;
    genomes4RunEAHelitron;
    genomes4MiteFinder;
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

    publishDir "${params.outdir}/noncoding/${name}"

    tag { name }

    input:
    set val(name), file(fasta) from genomes4RunTRNAScan

    output:
    set val(name),
        file("${name}_trnascan.txt"),
        file("${name}_trnascan_ss.txt") into tRNAScanResults

    file "${name}_trnascan_iso.txt"
    file "${name}_trnascan_stats.txt"
    file "${name}_trnascan.bed"
    file "${name}_trnascan.fasta"

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

    publishDir "${params.outdir}/noncoding/${name}"

    tag { name }

    input:
    set val(name), file("ts.txt"), file("ss.txt") from tRNAScanResults

    output:
    set val(name), file("${name}_trnascan.gff3") into tRNAScanGFF

    script:
    """
    gffpal trnascan2gff -o "${name}_trnascan.gff3" ts.txt ss.txt
    """
}


process pressRfam {

    label "infernal"
    label "small_task"

    when:
    params.structrnafinder

    input:
    file "Rfam.cm" from rfam

    output:
    file "out" into pressedRfam

    script:
    """
    mkdir out
    cp -L Rfam.cm out/Rfam.cm
    cmpress -F out/Rfam.cm
    """
}


// TODO replace this with regular infernal
process runStructRNAFinder {

    label "structrnafinder"
    label "big_task"

    publishDir "${params.outdir}/noncoding/${name}"

    tag { name }

    when:
    params.structrnafinder

    input:
    set val(name), file(fasta) from genomes4RunStructRNAFinder
    file "rfamdb" from pressedRfam

    output:
    set val(name),
        file("${name}_structrnafinder.tsv") into structRNAfinderResults
    file "${name}_structrnafinder.txt"
    file "html"
    file "img"

    script:
    """
    # Do something about truncating fasta headers
    structRNAfinder \
      -i "${fasta}" \
      -d rfamdb/Rfam.cm \
      -r \
      -c ${task.cpus} \
      --method cmsearch \
      --tblout "${name}_structrnafinder.tsv" \
      --output "${name}_structrnafinder.txt"
    """
}


/*
 * RNAmmer
 * doi: 10.1093/nar/gkm160
 * url:
 */
process runRNAmmer {

    label "rnammer"
    label "small_task"

    publishDir "${params.outdir}/noncoding/${name}"

    tag { name }

    when:
    params.rnammer

    input:
    set val(name), file(fasta) from genomes4RunRNAmmer

    output:
    set val(name), file("${name}_rnammer.gff2") into rnammerResults
    file "${name}_rnammer.hmmreport"

    """
    rnammer \
      -S euk \
      -m lsu,ssu,tsu \
      -gff "${name}_rnammer.gff2" \
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

    publishDir "${params.outdir}/noncoding/${name}"

    tag { name }

    input:
    set val(name), file("rnammer.gff2") from rnammerResults

    output:
    set val(name), file("${name}_rnammer.gff3") into rnammerGFF

    script:
    """
    gffpal rnammer2gff -k euk -o "${name}_rnammer.gff3" rnammer.gff2
    """
}


/*
 * OcculterCut
 * doi: 10.1093/gbe/evw121
 */
process runOcculterCut {

    label "occultercut"
    label "small_task"
    tag { name }
    publishDir "${params.outdir}/composition/${name}"

    input:
    set val(name), file(fasta) from genomes4RunOcculterCut

    output:
    file "${name}_occultercut_regions.gff3"
    file "${name}_occultercut.png"
    file "${name}_occultercut_composition_gc.txt"
    file "${name}_occultercut_my_genome.txt"
    file "${name}_occultercut_grouped_regions.gff3" optional true
    file "${name}_occultercut_nuc_frequencies.R*" optional true

    script:
    """
    OcculterCut -f "${fasta}"


    sed -i '1i set terminal pngcairo size 900,600 enhanced font "Helvetica,20"' plot.plt
    sed -i '1a set output "plot.png"' plot.plt
    gnuplot plot.plt

    mv plot.png "${name}_occultercut.png"

    mv compositionGC.txt "${name}_occultercut_composition_gc.txt"
    mv regions.gff3 "${name}_occultercut_regions.gff3"
    mv myGenome.txt "${name}_occultercut_my_genome.txt"

    if [ -e groupedRegions.gff3 ]
    then
      mv groupedRegions.gff3 "${name}_occultercut_grouped_regions.gff3"
    fi

    for f in nuc_frequencies.R*
    do
      mv "\${f}" "${name}_occultercut_\${f}"
    done
    """
}


//
// Find repeats/TEs
//


/*
 * RepeatModeller
 * url: http://www.repeatmasker.org/RepeatModeler/
 */
process runRepeatModeller {

    label "repeatmasker"
    label "medium_task"
    tag "${name}"
    publishDir "${params.outdir}/repeatmodeller"

    input:
    set val(name), file(fasta) from genomes4RunRepeatModeller
    file "rmlib" from repbase

    output:
    file "${name}_repeatmodeller_consensus.fasta"
    set val(name),
        file("${name}_repeatmodeller_msa.stk") into repeatModellerResults

    script:
    """
    export RM_LIB="\${PWD}/rmlib"

    # Where will this be placed?
    BuildDatabase -name "${name}" -engine ncbi "${fasta}"

    # Can we split this over scaffolds?
    RepeatModeler -engine ncbi -pa ${task.cpus} -database ${name} >& run.out

    mv "${name}-families.fa" "${name}_repeatmodeller_consensus.fasta"
    mv "${name}-families.stk" "${name}_repeatmodeller_msa.stk"
    """
}


/*
 * RepeatMasker
 * url: http://www.repeatmasker.org/RMDownload.html
 */
process runRepeatMaskerRepbase {

    label "repeatmasker"
    label "medium_task"
    tag "${name}"
    publishDir "${params.outdir}/repeatmasker"

    input:
    set val(name), file(fasta) from genomes4RunRepeatMaskerRepbase
    file "rmlib" from repbase

    output:
    set val(name), file("${name}") into repeatMasked

    script:
    """
    export RM_LIB="\${PWD}/rmlib"

    RepeatMasker \
      -e ncbi \
      -species "${params.rmspecies}" \
      -pa "${task.cpus}" \
      -xsmall \
      -dir "${name}" \
      "${fasta}"
    """
}


/*
 * LTRHarvest
 * doi: 10.1186/1471-2105-9-18
 */
process runLtrHarvest {

    label "genometools"
    label "small_task"

    tag "${name}"

    input:
    set val(name), file(fasta) from genomes4RunLtrHarvest

    output:
    set val(name), file(fasta), file("${fasta}*"), file("${name}_ltrharvest.gff3") into ltrHarvestResults
    set val(name), file("${name}_ltrharvest.fasta")
    set val(name), file("${name}_ltrharvest_inner.fasta")

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

    gt ltrharvest \
      -seqids yes \
      -index "${fasta}" \
      -gff3 "${name}_ltrharvest_tmp.gff3" \
      -out "${name}_ltrharvest.fasta" \
      -outinner "${name}_ltrharvest_inner.fasta"

    gt gff3 \
      -tidy \
      -sort \
      -retainids \
      "${name}_ltrharvest_tmp.gff3" \
    > "${name}_ltrharvest.gff3"
    """
}


/*
 * HMMs from gypsydb are in HMMER2 format, we just convert all to the
 * version of hmmer we're using.
 */
process fixLtrHMMs {

    label "hmmer3"
    label "small_task"

    input:
    file "in/*" from ltrHMMs.collect()

    output:
    file "out/*" into fixedLtrHMMs

    script:
    """
    mkdir -p out

    for f in in/*;
    do
        hmmconvert "\${f}" > "out/\$(basename \${f})"
    done
    """
}


/*
 * LTRDigest
 * doi: 10.1093/nar/gkp759
 *
 * TODO: filter out incomplete LTRs?
 * This should remove many false positives, but might might exclude things.
 */
process runLtrDigest {

    label "genometools"
    label "small_task"

    tag "${name}"

    input:
    set val(name),
        file(fasta),
        file("*"),
        file("${name}_ltrharvest.gff3"),
        file("hmms/*hmm") from ltrHarvestResults
            .combine(fixedLtrHMMs.toList())

    script:
    """
    mkdir -p tmp
    TMPDIR="\${PWD}/tmp" \
    gt -j "${task.cpus}" ltrdigest \
      -hmms hmms/*hmm \
      -outfileprefix "${name}_ltrdigest" \
      -matchdescstart \
      -seqfile "${fasta}" \
      "${name}_ltrharvest.gff3" \
    > ltrdigest.gff3

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
    tag "${name}"

    input:
    set val(name), file(fasta) from genomes4RunEAHelitron

    output:
    file "${name}.3.txt"
    file "${name}.5.txt"
    file "${name}.5.fa"
    file "${name}.u*.fa"
    file "${name}.d*.fa"
    file "${name}.gff3"
    file "${name}.bed"
    file "${name}.len.txt"

    script:
    three_prime_fuzzy_level = 4
    upstream_length = 3000 // default
    downstream_length = 500 // default

    """
    # -p "${task.cpus}"

    EAHelitron \
      -o "${name}" \
      -r "${three_prime_fuzzy_level}" \
      -u "${upstream_length}" \
      -d "${downstream_length}" \
      "${fasta}"
    """
}


/*
 * HelitronScanner
 * doi: 10.1073/pnas.1410068111
 * url: https://sourceforge.net/projects/helitronscanner/
 *
 * Todo, parse outputs into gff or bed-like file.
 * look at adding new regular expressions to lcvs files.
 * Figure out how to specify path of .jar archive
 * (probably set environment variable?)
process runHelitronScanner {

    label "helitronscanner"
    label "small_task"

    input:
    set val(name), file(genome) from genomes4HelitronScanner
    set file("heads.lcvs"), file("tails.lcvs") from

    output:
    set val(name), file("${genome.baseName}_fwd_pairs.txt"),
        file("${genome.baseName}_rev_pairs.txt") into helitronScannerLocations
    set val(name), file("${genome.baseName}_helitrons.fasta") into helitronScannerSeqs

    script:
    threshold = 5
    min_size = 200
    max_size = 50000
    """
    java -jar HelitronScanner.jar scanHead \
      -threads_LCV ${task.cpus} \
      -buffer_size 0 \
      -lcv_filepath heads.lcvs \
      -genome ${genome} \
      -output fwd_head_positions.txt \
      -overlap 50 \
      -threshold 1

    java -jar HelitronScanner.jar scanTail \
      -threads_LCV ${task.cpus} \
      -buffer_size 0 \
      -lcv_filepath tails.lcvs \
      -genome ${genome} \
      -output fwd_tail_positions.txt \
      -overlap 50 \
      -threshold 1

    java -jar HelitronScanner.jar pairends \
      -head_score fwd_head_positions.txt \
      -tail_score fwd_tail_positions.txt \
      -output "${genome.baseName}_fwd_pairs.txt" \
      -head_threshold ${threshold} \
      -tail_threshold ${threshold} \
      -helitron_len_range ${min_size}:${max_size}

    java -jar HelitronScanner.jar draw \
      -pscore "${genome.baseName}_fwd_pairs.txt" \
      -genome ${genome} \
      -output "${genome.baseName}_fwd" \
      --pure

    # Get rc matches
    java -jar HelitronScanner.jar scanHead \
      -threads_LCV ${task.cpus} \
      -buffer_size 0 \
      -lcv_filepath heads.lcvs \
      -genome ${genome} \
      -output rev_head_positions.txt \
      -overlap 50 \
      -threshold 1 \
      --rc_mode

    java -jar HelitronScanner.jar scanTail \
      -threads_LCV ${task.cpus} \
      -buffer_size 0 \
      -lcv_filepath tails.lcvs \
      -genome ${genome} \
      -output rev_tail_positions.txt \
      -overlap 50 \
      -threshold 1 \
      --rc_mode

    java -jar HelitronScanner.jar pairends \
      -head_score rev_head_positions.txt \
      -tail_score rev_tail_positions.txt \
      -output ${genome.baseName}_rev_pairs.txt \
      -head_threshold ${threshold} \
      -tail_threshold ${threshold} \
      -helitron_len_range ${min_size}:${max_size} \
      --rc_mode

    java -jar HelitronScanner.jar draw \
      -pscore ${genome.baseName}_rev_pairs.txt \
      -genome ${genome} \
      -output "${genome.baseName}_rev" \
      --pure

    cat "${genome.baseName}_fwd.hel.fa" "${genome.baseName}_rev.hel.fa" \
      > ${genome.baseName}_helitrons.fasta
    """
}
 */


/*
 * MITEfinderII
 * doi: 10.1186/s12920-018-0418-y
 * url: https://github.com/screamer/miteFinder
 *
 * Compare performance with MiteTracker if can resolve the lousy install process.
 */
process runMiteFinder {
    label "mitefinder"
    tag { name }

    input:
    set val(name), file(genome) from genomes4MiteFinder
    file "pattern_scoring.txt" from miteFinderProfiles

    output:
    set val(name), file("${name}_mitefinder.fasta") into miteFinderFasta

    script:
    threshold = 0.5

    """
    miteFinder \
      -input "${genome}" \
      -output "${name}_mitefinder.fasta" \
      -pattern_scoring pattern_scoring.txt \
      -threshold "${threshold}"
    """
}


/*
 * Combine and remove redundancy.
 * Meshclust looks like good tool for clustering, (alternative to cdhit-cdna)
 * Possibly need to look at overlapping genomic annotations?
 * Filtering out false positives? negative Database matches?
 */


/*
 * Multiple sequence alignment and annotation
 * Mafft is always solid aligner choice.
 * Paste probably good option for classification.
 * Classification might go hand in hand with filtering false positives?
 * Maybe fp is first step is to exclude members of clusters,
 *  whereas here it is to exclude entire clusters?
 */
