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
    tuple val(sample_id), path(reads) from ch_18S_input

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
    tuple val(sample_id), path("${sample_id}.fasta"), path("${sample_id}.count"), path("${sample_id}.summary.logfile") into ch_merge_pcr_reads_18S_out
    tuple val(sample_id), path("${sample_id}.fasta"), path("${sample_id}.count") into ch_make_mmseqs_query_dbs

    script:
    """
    cat $fasta_fwd $fasta_rev > ${sample_id}.fasta
    cat $names_fwd $names_rev > ${sample_id}.names
    mothur "#count.seqs(name=${sample_id}.names, processors=${task.cpus})"
    rm mothur.*.logfile
    mothur "#summary.seqs(fasta=${sample_id}.fasta, name=${sample_id}.names, processors=${task.cpus})"
    mv mothur.*.logfile ${sample_id}.summary.logfile
    mv ${sample_id}.count_table ${sample_id}.count
    """
}

// Create mmseqs db of the fasta files
// This is very memory hungry. Memory consumption is estimated at 75G per search, but the below parameters
// (memory 80G and cpus 20) seems to be working well
// htop shows 80% CPU and 615G memory.
process make_mmseqs_tax_table{
    tag {sample_id}
    container "soedinglab/mmseqs2:latest"
    cpus 20
    memory "80G"
    if (workflow.containerEngine == 'docker'){
        containerOptions "-v ${params.tmp_parent_dir}:${params.tmp_parent_dir}"
    }
    publishDir "taxonomy_tables/${sample_id}/", mode: "copy"

    input:
    tuple val(sample_id), path(fasta), path(count), path(silva_db_base), path(silva_db_indices) from ch_make_mmseqs_query_dbs.combine(Channel.fromPath(params.nt_18S_path)).combine(Channel.fromPath("${params.nt_18S_path}{.,_}*").collect().map{[it]})

    output:
    tuple val(sample_id), path(count), path("${sample_id}.taxonomyResult.tsv") into ch_make_krona_input

    script:
    """
    # Create the mmseqs query db
    mmseqs createdb ${fasta} ${sample_id}.queryDB
    # Remove the temporary directory if it already exists
    rm -r ${params.tmp_parent_dir}/${sample_id} || true
    # Do the taxonomy search
    mmseqs taxonomy --threads ${task.cpus} --search-type 3 --tax-lineage 1 ${sample_id}.queryDB $silva_db_base ${sample_id}.taxonomyResult ${params.tmp_parent_dir}/${sample_id}
    # Clean Up: Remove the temp dir
    rm -r ${params.tmp_parent_dir}/${sample_id}
    # Create the taxonomy tsv
    mmseqs createtsv ${sample_id}.queryDB ${sample_id}.taxonomyResult ${sample_id}.taxonomyResult.tsv
    """
}

// The awk that ships with many of the docker images is mawk
// and the for loops are not working properly so I have created a
// very basic "gawk" image to do the awk commands in.
// For any sequeneces that weren't mapped, we want to summate them and add
// them to the end of the krona input file.
process make_krona_input{
    tag {sample_id}
    container "didillysquat/gawk:v5.0.1-ubuntu_20210225"
    if (workflow.containerEngine == 'docker'){
            containerOptions '-u $(id -u):$(id -g)'
        }
    cpus 1
    publishDir "krona_input/${sample_id}/", mode: "copy"
    
    input:
    tuple val(sample_id), path(count_table), path(tax_table) from ch_make_krona_input

    output:
    tuple val(sample_id), path("${sample_id}.krona.input.txt") into ch_make_krona_plot,ch_make_summary_krona_inputs
    
    shell:
    '''
    # Use the taxonomy tsv and the count file to create input for Krona figure
    awk 'BEGIN{FS="\t"} NR==FNR {if(NR>1){abund_map[$1]=$2; next}} {new_tax=""; split($5, tax_a, ";"); for(i=1; i<=length(tax_a); i++){if(substr(tax_a[i],1,1)!="-"){new_tax=new_tax "\t" tax_a[i]};}; if(new_tax!=""){tax_map[new_tax]+=abund_map[$1]};} END {for(tax in tax_map){split_tax = tax; gsub(";","\t", split_tax); print tax_map[tax] "\t" split_tax;}}' !{count_table} !{tax_table} > !{sample_id}.krona.input.txt
    # Check to see if there are non-matched sequences and summate these and add to end of input if so
    awk 'BEGIN{FS="\t"; i=0; unmatched=0} NR==FNR {matched_array[$1]=i; i++; next;}{if(FNR>1){if(!($1 in matched_array)){unmatched+=$2}}} END {print unmatched}' !{tax_table} !{count_table} >> !{sample_id}.krona.input.txt
    '''
}

process make_krona_plot{
    tag {sample_id}
    container "biocontainers/krona:v2.7.1_cv1"
    if (workflow.containerEngine == 'docker'){
            containerOptions '-u $(id -u):$(id -g)'
        }
    cpus 1
    publishDir "krona_output/${sample_id}/", mode: "copy"
    
    input:
    tuple val(sample_id), path(krona_input) from ch_make_krona_plot

    output:
    tuple val(sample_id), path("${sample_id}.krona.html") into ch_make_krona_out
    
    script:
    """
    # Make the Krona figure using the input
    ktImportText ${krona_input}
    mv text.krona.html ${sample_id}.krona.html
    """    
}

// Make kronas that uses per sample normalised abundances
// Make one that is just cnidaria broken down. For this one we only keep the most abundant p_Cnidaria sequence per sample. We give this a value of 1.
// As such the output krona plot esentially gives a split of the taxa as a sample count.
// Make one that is split into cnidaria, symbiodiniaceae and other. For this one we have done a per sample normalisation and added up these realtive abundances
// As such the output gives a total sequence abundance of the given taxa groups across the whole dataset.
process make_multi_sample_host_krona{
    cpus 1
    publishDir "krona_summaries/", mode: "copy"

    input:
    path(krona_input) from ch_make_summary_krona_inputs.collect{it[1]}

    output:
    tuple path("cnidaria.krona.input.txt"), path("categories.krona.input.txt") into ch_make_summary_krona_plots

    shell:
    '''
    # Because krona automatically summates lines that are the same
    # We can simply cat together all of the krona inputs once
    # they have been per sample normalised
    # TODO per sample, do normalisation, filter into cnidarian and three categories
    for KINPUT in *krona.input.txt
    do
    # Pull out the most abundant cnidarian sequence and give it value of 1
    awk 'BEGIN{FS="\t";big=0;} NR==FNR {seq_tot+=$1; next;} {if(match($0, /p_Cnidaria/)){if($1/seq_tot > big){big=$1/seq_tot; new_string=1; for(i=2; i<=NF; i++){new_string = new_string "\t" $i}};}} END{print new_string >> "cnidaria.krona.input.txt"}' $KINPUT $KINPUT
    # Categorize into Cnidaria, Symbiodiniaceae, Other
    awk 'BEGIN{FS="\t"} NR==FNR {seq_tot+=$1; next;} {if(match($0, /p_Cnidaria/)){new_string=$1/seq_tot "\tp_Cnidaria"}else if(match($0, /f_Symbiodiniaceae/)){new_string=$1/seq_tot "\tf_Symbiodiniaceae"}else{new_string=$1/seq_tot "\tOther"}; print new_string >> "categories.krona.input.txt";}' $KINPUT $KINPUT
    done
    '''
}

process make_krona_summary_plots{
    tag {sample_id}
    container "biocontainers/krona:v2.7.1_cv1"
    if (workflow.containerEngine == 'docker'){
            containerOptions '-u $(id -u):$(id -g)'
        }
    cpus 1
    publishDir "krona_summaries", mode: "copy"
    
    input:
    tuple path(cnidarian_input), path(category_input) from ch_make_summary_krona_plots

    output:
    tuple path("cnidaria.krona.html"), path("categories.krona.html") into ch_make_krona_summary_out
    
    script:
    """
    # Make the Krona figure using the input
    ktImportText ${cnidarian_input}
    mv text.krona.html cnidaria.krona.html
    ktImportText ${category_input}
    mv text.krona.html categories.krona.html
    """    
}