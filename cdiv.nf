#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

params.input_seq_dir = "/home/humebc/projects/tara/cdiv/inputs/seq_download"
params.tn_map_csv_path = "/home/humebc/projects/tara/cdiv/inputs/tn_map.csv"

// Setup up the three input channels
// This will create a channel for each marker that is a list of tuples where the first
// item is the sample-id and the second is a tuple of the R1 and R2 read in that order
(channel_list_16S, channel_list_18S, channel_list_ITS) = new File(params.tn_map_csv_path).withReader { reader ->
    def tn_map_csv = parseCsv(reader)
    def tup_list_16S = []
    def tup_list_18S = []
    def tup_list_ITS = []
    for(line in tn_map_csv) {
        if (line.primers.contains("16S")){
            tup_list_16S << ["TARA_${line.barcode}", [file([params.input_seq_dir, line.R1.tokenize("/").last()].join(File.separator)), file([params.input_seq_dir, line.R2.tokenize("/").last()].join(File.separator))]]
        }else if(line.primers.contains("18S")){
            tup_list_18S << ["TARA_${line.barcode}", [file([params.input_seq_dir, line.R1.tokenize("/").last()].join(File.separator)), file([params.input_seq_dir, line.R2.tokenize("/").last()].join(File.separator))]]
        }else if(line.primers.contains("ITS")){
            tup_list_ITS << ["TARA_${line.barcode}", [file([params.input_seq_dir, line.R1.tokenize("/").last()].join(File.separator)),file([params.input_seq_dir, line.R2.tokenize("/").last()].join(File.separator))]]
        }
    }
    return [tup_list_16S, tup_list_18S, tup_list_ITS]
}

ch_16S_input = Channel.fromList(channel_list_16S)
ch_18S_input = Channel.fromList(channel_list_18S)
ch_ITS_input = Channel.fromList(channel_list_ITS)

process make_oligos_file{
    tag "make_oligos"

    input:
    tuple val(fwd_primer), val(rev_primer) from Channel.of(["TTGTACACACCGCCC", "CCTTCYGCAGGTTCACCTAC"])
    tuple val(fwd_primer_rc), val(rev_primer_rc) from Channel.of(["CCTTCYGCAGGTTCACCTAC", "TTGTACACACCGCCC"])

    output:
    path("18S.oligos.oligos") into ch_18S_oligos
    path("18S.oligos.contigs.oligos") into ch_make_contigs_oligos
    path("18S.oligos.rc.oligos") into ch_18S_rc_oligos

    script:
    """
    printf "forward\t${fwd_primer}\nreverse\t${rev_primer}" > 18S.oligos.oligos
    printf "primer\t${fwd_primer}\t${rev_primer}" > 18S.oligos.contigs.oligos
    printf "forward\t${fwd_primer_rc}\nreverse\t${rev_primer_rc}" > 18S.oligos.rc.oligos
    """
}

// // The currently available biocontainers/mothur image is broken
// // (make.contigs doesn't work with gzipped fastq files)
// // Working with conda instead
// process make_contigs_18S{
//     tag {sample_id}
//     container "didillysquat/mothur:v1.45.0-ubuntu_20210225"
//     cpus 1

//     input:
//     tuple val(sample_id), path(reads), path(oligos) from ch_18S_input.take(4).combine(ch_make_contigs_oligos)

//     output:
//     tuple val(sample_id), path("${sample_id}.contigs.fasta") into ch_make_contigs_unique
    
//     script:
//     """
//     mothur "#make.contigs(ffastq=${reads[0]}, rfastq=${reads[1]}, processors=${task.cpus}, oligos=${oligos}, pdiffs=4)"
//     mv *fastq.trim.contigs.fasta ${sample_id}.contigs.fasta
//     """
// }

process make_contigs_18S{
    tag {sample_id}
    container "didillysquat/mothur:v1.45.0-ubuntu_20210225"
    cpus 1

    input:
    tuple val(sample_id), path(reads) from ch_18S_input.take(4)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.fasta") into ch_make_contigs_unique
    
    script:
    """
    mothur "#make.contigs(ffastq=${reads[0]}, rfastq=${reads[1]}, processors=${task.cpus})"
    mv *fastq.trim.contigs.fasta ${sample_id}.contigs.fasta
    """
}

process unique_contigs_18S{
    tag {sample_id}
    container "didillysquat/mothur:v1.45.0-ubuntu_20210225"
    cpus 1

    input:
    tuple val(sample_id), path(contig_fastq) from ch_make_contigs_unique

    output:
    tuple val(sample_id), path("${sample_id}.contigs.unique.fasta"), path("${sample_id}.contigs.names") into ch_contigs_summary_18S,ch_pcr_18S

    script:
    """
    mothur "#unique.seqs(fasta=${contig_fastq}, processors=${task.cpus})"
    """
}

process unique_contigs_summary_18S{
    tag {sample_id}
    container "didillysquat/mothur:v1.45.0-ubuntu_20210225"
    cpus 1
    publishDir "unique_contig_summary", mode: "copy"

    input:
    tuple val(sample_id), path(fasta), path(names) from ch_contigs_summary_18S

    output:
    tuple val(sample_id), path("${sample_id}.contigs.unique.summary"), path("${sample_id}.contigs.unique.summary.logfile") into ch_contigs_summary_out_18S

    script:
    """
    mothur "#summary.seqs(fasta=${fasta}, name=${names}, processors=${task.cpus})"
    mv *.logfile ${sample_id}.contigs.unique.summary.logfile
    """
}

process pcr_18S_fwd{
    tag {sample_id}
    container "didillysquat/mothur:v1.45.0-ubuntu_20210225"
    cpus 1
    publishDir "fwd_pcr_out", mode: "copy"

    input:
    tuple val(sample_id), path(fasta), path(names), path(oligos) from ch_pcr_18S.combine(ch_18S_oligos)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.unique.pcr.unique.fasta"), path("${sample_id}.contigs.unique.pcr.names") into ch_pcr_summary_18S_fwd,ch_merge_pcrs_18S_fwd
    tuple val(sample_id), path("${sample_id}.contigs.unique.scrap.pcr.fasta"), path(names) into ch_pcr_18S_rev

    script:
    """
    mothur "#pcr.seqs(fasta=${fasta}, name=${names}, oligos=${oligos}, processors=${task.cpus}, pdiffs=1, rdiffs=1)"
    mothur "#unique.seqs(fasta=${sample_id}.contigs.unique.pcr.fasta, name=${sample_id}.contigs.pcr.names)"
    """
}

// Here we are performing a pcr on the scrap sequences to catch those sequences that were reverse complement
// We do this with the olio file reversed. We then do a reverse complement to get all seqs in the same orientation
// Finally, we perform a unique.seqs.
// We ivnested a lot of time getting the awk statement to compile properly with nextflow.
// The key was to escape the dollar signs AND the back slash. We were having issues getting the shell format to work.
process pcr_18S_rev{
    tag {sample_id}
    container "didillysquat/mothur:v1.45.0-ubuntu_20210225"
    cpus 1
    publishDir "rev_pcr_out", mode: "copy"

    input:
    tuple val(sample_id), path(fasta), path(names), path(oligos) from ch_pcr_18S_rev.combine(ch_18S_rc_oligos)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.unique.scrap.pcr.clean.pcr.unique.fasta"), path("${sample_id}.contigs.unique.scrap.pcr.clean.pcr.names") into ch_pcr_summary_18S_rev,ch_merge_pcr_18S_rev

    script:
    """
    awk '{if(match(\$1, /^>/)){gsub(/\\|[frt]*/, "", \$1); print \$1;}else{print \$0}}' ${fasta} > ${sample_id}.contigs.unique.scrap.pcr.clean.fasta
    mothur "#pcr.seqs(fasta=${sample_id}.contigs.unique.scrap.pcr.clean.fasta, name=${names}, oligos=${oligos}, processors=${task.cpus}, pdiffs=1, rdiffs=1)"
    mv ${sample_id}.contigs.pcr.names ${sample_id}.contigs.pcr.rc.names
    mothur "#reverse.seqs(fasta=${sample_id}.contigs.unique.scrap.pcr.clean.pcr.fasta)"
    mv ${sample_id}.contigs.unique.scrap.pcr.clean.pcr.rc.fasta ${sample_id}.contigs.unique.scrap.pcr.clean.pcr.fasta
    mothur "#unique.seqs(fasta=${sample_id}.contigs.unique.scrap.pcr.clean.pcr.fasta, name=${sample_id}.contigs.pcr.rc.names)"
    """
}

process pcr_summary_18S{
    tag {sample_id}
    container "didillysquat/mothur:v1.45.0-ubuntu_20210225"
    cpus 1
    publishDir "pcr_summary", mode: "copy"

    input:
    tuple val(sample_id), path(fasta_fwd), path(names_fwd), path(fasta_rev), path(names_rev) from ch_pcr_summary_18S_fwd.join(ch_pcr_summary_18S_rev)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.unique.pcr.unique.summary"), path("${sample_id}.contigs.unique.pcr.unique.summary.logfile") into ch_pcr_summary_18S_out_fwd
    tuple val(sample_id), path("${sample_id}.contigs.unique.scrap.pcr.clean.pcr.unique.summary"), path("${sample_id}.contigs.unique.scrap.pcr.clean.pcr.unique.summary.logfile") into ch_pcr_summary_18S_out_rev

    script:
    """
    mothur "#summary.seqs(fasta=${fasta_fwd}, name=${names_fwd}, processors=${task.cpus})"
    mv mothur.*.logfile ${sample_id}.contigs.unique.pcr.unique.summary.logfile

    mothur "#summary.seqs(fasta=${fasta_rev}, name=${names_rev}, processors=${task.cpus})"
    mv mothur.*.logfile ${sample_id}.contigs.unique.scrap.pcr.clean.pcr.unique.summary.logfile
    """
}

process merge_pcr_reads_18S{
    tag {sample_id}
    container "didillysquat/mothur:v1.45.0-ubuntu_20210225"
    cpus 1
    publishDir "post_pcr", mode: "copy"
    
    input:
    tuple val(sample_id), path(fasta_fwd), path(names_fwd), path(fasta_rev), path(names_rev) from ch_merge_pcrs_18S_fwd.join(ch_merge_pcr_18S_rev)

    output:
    tuple val(sample_id), path("${sample_id}.fasta"), path("${sample_id}.count"), path("${sample_id}.summary.logfile") into ch_merge_pcr_reads_18S

    script:
    """
    cat $fasta_fwd $fasta_rev > ${sample_id}.fasta
    cat $names_fwd $names_rev > ${sample_id}.names
    mothur "#count.seqs(name=${sample_id}.names)"
    rm mothur.*.logfile
    mothur "#summary.seqs(fasta=${sample_id}.fasta, name=${sample_id}.names)"
    mv mothur.*.logfile ${sample_id}.summary.logfile
    mv ${sample_id}.count_table ${sample_id}.count
    """



}