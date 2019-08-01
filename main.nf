#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # pante


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

params.structrnafinder = false
params.rfam = false
params.rfam_url = "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
params.rnammer = false


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


genomes.into {
    genomes4RunTRNAScan;
    genomes4RunStructRNAFinder;
    genomes4RunRNAmmer;
    genomes4RunOcculterCut;
}


//
// Finding non-coding RNA
//

process runTRNAScan {

    label "trnascan"
    label "medium_task"

    publishDir "${params.outdir}/noncoding/${name}"

    tag { name }

    input:
    set val(name), file(fasta), file(faidx) from genomes4RunTRNAScan

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
      --log trna.log \
      --thread "${task.cpus}" \
      "${fasta}"
    """
}


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
    trnascan2gff -o "${name}_trnascan.gff3" ts.txt ss.txt
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


process runStructRNAFinder {

    label "structrnafinder"
    label "big_task"

    publishDir "${params.outdir}/noncoding/${name}"

    tag { name }

    when:
    params.structrnafinder

    input:
    set val(name), file(fasta), file(faidx) from genomes4RunStructRNAFinder
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


process runRNAmmer {

    label "rnammer"
    label "small_task"

    publishDir "${params.outdir}/noncoding/${name}"

    tag { name }

    when:
    params.rnammer

    input:
    set val(name), file(fasta), file(faidx) from genomes4RunRNAmmer

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
    rnammer2gff -k euk -o "${name}_rnammer.gff3" rnammer.gff2
    """
}


process runOcculterCut {

    label "occultercut"
    label "small_task"

    publishDir "${params.outdir}/composition/${name}"

    tag { name }

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


// Run


process runRepeatMasker {

    label "repeatmasker"
    label "medium_task"

    publishDir "${params.outdir}/repeatmasker"

    tag "${name}"

    input:
    set val(name), file(genome) from genomes
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
      "${genome}"
    """
}

/*
 * De-novo repeat assemblers
 */

// REPdenovo

/*
 * De-novo repeat finders
 */


// Repeat modeller
/*
process runRepeatModeller {
    label "repeatmasker"

    input:
    set val(name), file(genome) from genomes4RepeatModeller


    """
    # Where will this be placed?
    BuildDatabase -name ${name} -engine ncbi ${genome}

    # Can we split this over scaffolds?
    RepeatModeler -engine ncbi -pa ${task.cpus} -database ${name} >& run.out
    """
}
*/

// CARP ? Might need to pipeline this myself, currently exists as separate tools and description of workflow.
// MGEScan-non-ltr

/*
 * Homology based repeat finders
 */

// Repeat masker

/*
 * LTR finders
 */


// LTR_harvest + digest // See ltr_retriever paper for "good" parameters
// LTR_retriever post-processor for ltr_harvest. Possibly need to exclude
// because of dependence on RepeatMasker? Might be tricky to get running.

// LTR_detector
// Doesn't seem to find many in stago?
/*
process runLtrDetector {
    label "ltrdetector"

    """
    # Split multifasta into directory with one chrom per file

    LtrDetector \
        -fasta fasta_dir/ \
        -destDir output_dir \
        -id 70 \
        -nThreads ${task.cpus} \
        -nested
    """
}
*/

/*
process processLtrDetectorResults {

    """
    Combine bed-like files.
    """
}
*/

/*
 * SINE finders
 */

// SINE_Scan


/*
 * MITE finders
 */

/*
 * MITEfinderII
 *
 * TODO: Factor out "pattern_scoring.txt" as input parameter.
 * Can't see a way to make default to environment variable in docker image.
 * Possibly redistribute with pipeline?
 */
/*
process runMiteFinder {
    label "mitefinder"
    tag { name }

    input:
    set val(name), file(genome) from genomes4MiteFinder

    output:
    set val(name), file("mf.fasta") into miteFinderFasta

    """
    miteFinder \
      -input ${genome} \
      -output mf.fasta \
      -pattern_scoring /opt/mitefinder/profile/pattern_scoring.txt \
      -threshold 0.5
    """
}
*/

/*
process processMiteFinder {
    label "python3"
    tag { name }

    // Convert fasta into bed or gff.
}
*/

// MITEHunter

/*
 * Helitron finders

// HelitronScanner
// Todo, parse outputs into gff or bed-like file.
// look at adding new regular expressions to lcvs files.
// Figure out how to specify path of .jar archive
// (probably set environment variable?)
process runHelitronScanner {
    label "helitronscanner"

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
