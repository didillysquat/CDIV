#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

params.input_seq_dir = "/home/humebc/projects/tara/cdiv/inputs/seq_download"
params.tn_map_csv_path = "/home/humebc/projects/tara/cdiv/inputs/tn_map.csv"

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
ch_18S_input = Channel.fromList(channel_list_16S)
ch_ITS_input = Channel.fromList(channel_list_16S)

ch_16S_input.view()