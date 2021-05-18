"""
Script to answer some key quesitons about the distribution of the 18S host sequences in the CDIV data
Question 1: From the 2400 samples, roughly how many unique most abundant host 18S sequences are there
"""

import os
import pandas as pd
import pickle
from collections import defaultdict
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'

class Questions:
    def __init__(self):
        self.tax_table_path = "/home/humebc/projects/tara/cdiv/output_backup_20210510/taxonomy_tables"
        self.post_pcr_path = "/home/humebc/projects/tara/cdiv/output_backup_20210510/post_pcr"
        self.resource_path = "/home/humebc/projects/tara/cdiv/script_resources"
        self.fig_dir = "/home/humebc/projects/tara/cdiv/figures"

        sample_list = []
        for path, directories, files in os.walk(self.tax_table_path):
            sample_list = directories
            break
        # Create a set of dictionaries for each sample
        if os.path.exists(os.path.join(self.resource_path, "sample_dict.p")):
            print("loading sample dict from pickle")
            sample_dict = pickle.load(open(os.path.join(self.resource_path, "sample_dict.p"), "rb"))
            print("done")
        else:
            # Dictionary of sample name as key with subdictionary as value. Subdict key is seq_name value is tuple of abundance_abs, abundans_rel, tax (as cnidarian, symbiodiniaceae, other), sequence
            # Corressponding dictionary where subdict is seq_name, value is tax as cnidarian, symbiodiniaceae or other
            # Corressponding dictionary where subdict is seq_name, value is nucleotide sequence
            # pickle these out
            print("no pickle of sample dict found. Creating from scratch.")
            sample_dict = {}
            for sample in sample_list:
                print(f"making subdict for {sample}")
                sub_dict = {}
                count_df = pd.read_table(os.path.join(self.tax_table_path, sample, f"{sample}.count"), index_col=0, sep="\t")
                count_dict = dict(count_df.iloc[:,0].items())
                count_tot = sum(count_dict.values())
                tax_df = pd.read_table(os.path.join(self.tax_table_path, sample, f"{sample}.taxonomyResult.tsv"), names=["seq_name", "tax_id", "tax_rank", "tax_val", "full_tax_str"], sep="\t")
                tax_df = tax_df.set_index("seq_name")
                tax_dict = {}
                tax_dict_full_str = {}
                for ind in tax_df.index:
                    tax_string = tax_df.at[ind, "full_tax_str"]
                    try:
                        if "p_Cnidaria" in tax_string:
                            tax_dict[ind] = "Cnidaria"
                        elif "f_Symbiodiniaceae" in tax_string:
                            tax_dict[ind] = "Symbiodiniaceae"
                        else:
                            tax_dict[ind] = "other"
                        tax_dict_full_str[ind] = tax_string
                    except TypeError:
                        tax_dict[ind] = "other"
                        tax_dict_full_str[ind] = "none"
                    
                with open(os.path.join(self.post_pcr_path, f"{sample}.fasta"), "r") as f:
                    fasta_list = [_.rstrip() for _ in f]
                seq_dict = {}
                for i, line in enumerate(fasta_list):
                    if i%2 == 0:
                        seq_dict[fasta_list[i].split("\t")[0][1:]] = fasta_list[i+1]
                # Here we have a dict for each of the sequence
                for sequence in count_df.index:
                    sub_dict[sequence] = (count_dict[sequence], count_dict[sequence]/count_tot, tax_dict[sequence], seq_dict[sequence], tax_dict_full_str[sequence])
                sample_dict[sample] = sub_dict
            foo = "bar"
            pickle.dump(sample_dict, open(os.path.join(self.resource_path, "sample_dict.p"), "wb"))
            
        # Now we can work on the question
        # First question, how many unique sequences make up the most abudant sequence
        print("\n\n\nCounting maj seqs")
        maj_seq_dd_path = os.path.join(self.resource_path, "maj_seq_dd.p")
        sample_to_maj_seq_dict_path = os.path.join(self.resource_path, "sample_to_maj_seq_dict.p")
        sample_to_maj_seq_name_dict_path = os.path.join(self.resource_path, "sample_to_maj_seq_name_dict.p")
        sequence_to_tax_string_dict_path = os.path.join(self.resource_path, "sequence_to_tax_string_dict.p")
        if os.path.exists(maj_seq_dd_path) and os.path.exists(sample_to_maj_seq_dict_path) and os.path.exists(sample_to_maj_seq_name_dict_path) and os.path.exists(sequence_to_tax_string_dict_path):
            maj_seq_dd = pickle.load(open(maj_seq_dd_path, "rb"))
            sample_to_maj_seq_dict = pickle.load(open(sample_to_maj_seq_dict_path, "rb"))
            sample_to_maj_seq_name_dict = pickle.load(open(sample_to_maj_seq_name_dict_path, "rb"))
            sequence_to_tax_string_dict = pickle.load(open(sequence_to_tax_string_dict_path, "rb"))
        else:
            maj_seq_dd = defaultdict(int)
            sample_to_maj_seq_dict = {}
            sample_to_maj_seq_name_dict = {}
            sequence_to_tax_string_dict = {}
            for sample, sub_dict in sample_dict.items():
                print(f"Getting maj seq for sample {sample}")
                maj_seq = [_[1][3] for _ in list(sorted(sub_dict.items(), key=lambda x: x[1][0], reverse=True)) if sub_dict[_[0]][2] == "Cnidaria"][0]
                maj_seq_name = [_[0] for _ in list(sorted(sub_dict.items(), key=lambda x: x[1][0], reverse=True)) if sub_dict[_[0]][2] == "Cnidaria"][0]
                maj_seq_tax_string = [_[1][4] for _ in list(sorted(sub_dict.items(), key=lambda x: x[1][0], reverse=True)) if sub_dict[_[0]][2] == "Cnidaria"][0]
                maj_seq_dd[maj_seq] += 1
                sample_to_maj_seq_dict[sample] = maj_seq
                sample_to_maj_seq_name_dict[sample] = maj_seq_name
                sequence_to_tax_string_dict[maj_seq] = maj_seq_tax_string
            pickle.dump(maj_seq_dd, open(maj_seq_dd_path, "wb"))
            pickle.dump(sample_to_maj_seq_dict, open(sample_to_maj_seq_dict_path, "wb"))
            pickle.dump(sample_to_maj_seq_name_dict, open(sample_to_maj_seq_name_dict_path, "wb"))
            pickle.dump(sequence_to_tax_string_dict, open(sequence_to_tax_string_dict_path, "wb"))

        # plot cululative proportion of the samples covered by the sequences
        sorted_maj_seqs = sorted(maj_seq_dd.items(), key = lambda x: x[1], reverse=True)
        cumulative_props = [0]
        cumulative_tot = 0
        tot = sum([_[1] for _ in sorted_maj_seqs])
        for seq, num_samples in sorted_maj_seqs:
            cumulative_tot += num_samples/tot
            cumulative_props.append(cumulative_tot)

        print(cumulative_props)     
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(range(len(cumulative_props)), cumulative_props, 'b--')
        ax.set_ylabel("Cumulative proportion of samples represented")
        ax.set_xlabel("Number of 18S sequences")
        plt.savefig(os.path.join(self.fig_dir, f"18S.cumul.prop.svg"))
        plt.savefig(os.path.join(self.fig_dir, f"18S.cumul.prop.png"), dpi=600)
        plt.close()

        # Next we want to make a tree or possibly a haplotype for visualisation.
        # To make a haplotype in R we will out put a df that has the sample as index,
        # the maj seq, and the taxonomy string of that seq
        haplo_df = pd.DataFrame(index=sample_list)
        haplo_df["maj_seq"] = [sample_to_maj_seq_dict[_] for _ in sample_list]
        haplo_df["tax_str"] = [sample_dict[_][sample_to_maj_seq_name_dict[_]][4] for _ in sample_list]
        haplo_df.to_csv(os.path.join(self.resource_path, "18S.haplotype.df.csv"))

        # We will also want to make a phylogenetic tree based on the maj 18S sequences
        # We will output a fasta that has the abundance of each of the samples built into its name
        fasta_list = []
        seq_count = 0
        tree_df_data = [[],[],[],[]]
        for seq, abund in sorted_maj_seqs:
            fasta_list.append(f">seq{seq_count}_{abund}")
            fasta_list.append(f"{seq}")
            tree_df_data[0].append(f"seq{seq_count}_{abund}")
            tree_df_data[1].append(seq)
            tree_df_data[2].append(abund)
            tree_df_data[3].append(sequence_to_tax_string_dict[seq])
            seq_count += 1
        with open(os.path.join(self.resource_path, "18S.maj.seqs.fasta"), "w") as f:
            for line in fasta_list:
                f.write(f"{line}\n")
        meta_df = pd.DataFrame(list(zip(*tree_df_data[1:])), columns=["seq", "abund", "tax_string"], index=tree_df_data[0])
        meta_df.to_csv(os.path.join(self.resource_path, "18S.maj.seqs.meta.df.csv"))
        # To assist in visualising the tree we will also output a set of meta data

        
        # We'd like to what the average proportion the most abundant sequence represents in a sample
        # we can plot this up as an ordered scatter
        maj_rel_abund_list = []
        cnidarian_rel_abund_list = []
        cnidarian_num_seqs_list = []
        for sample in sample_list:
            maj_host_seq_name = sample_to_maj_seq_name_dict[sample]
            rel_abund_of_host_maj_seq = sample_dict[sample][maj_host_seq_name][1]
            maj_rel_abund_list.append(rel_abund_of_host_maj_seq)
            # Also want to get a corresponding value of the relative abundance of all cnidarian sequences in the sample
            subdict = sample_dict[sample]
            cnidarian_abunds = [tup[1] for tup in subdict.values() if tup[2] == 'Cnidaria']
            cnidarian_rel_abund_list.append(sum(cnidarian_abunds))
            cnidarian_num_seqs_list.append(len(cnidarian_abunds))
        sorted_maj_rel_abund_list = sorted(maj_rel_abund_list, reverse=True)

        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.scatter(range(len(sorted_maj_rel_abund_list)), sorted_maj_rel_abund_list)
        ax.set_ylabel("Within sample proportion of most abundant host 18S sequence")
        ax.set_xlabel("Sample number")
        plt.savefig(os.path.join(self.fig_dir, f"18S.maj.seq.props.svg"))
        plt.savefig(os.path.join(self.fig_dir, f"18S.maj.seq.props.png"), dpi=600)
        plt.close()
        foo = "bar"

        # The above plot is a little worrying as it shows us that there are some samples that have the most abundant 18S as a low
        # abundance. The highest was about 0.5, but there were a few close and below the 0.1 mark
        # Next it would be good to see if these low abundance are due to there being a low level of Cnidarian sequnces in general
        # or whther they may be due to for example, a high richness and even diversity of 18S sequences (i.e. lots of lower level seqs)
        # To assess this we can plot up a scatter plot of Cnidarian proportion vs maj seq proportin
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.scatter(maj_rel_abund_list, cnidarian_rel_abund_list)
        ax.set_ylabel("Cnidarian rel abund")
        ax.set_xlabel("Most abund seq rel abund")
        plt.savefig(os.path.join(self.fig_dir, f"18S.cnid.vs.maj.props.svg"))
        plt.savefig(os.path.join(self.fig_dir, f"18S.cnid.vs.maj.props.png"), dpi=600)
        plt.close()

        # Then it is also useful to look at the number of unique 18S sequences as a product of the maj abund
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.scatter(maj_rel_abund_list, cnidarian_num_seqs_list, s=1)
        ax.set_ylabel("num cnidarian seqs")
        ax.set_xlabel("Most abund seq rel abund")
        plt.savefig(os.path.join(self.fig_dir, f"18S.cnid.vs.num.seqs.svg"))
        plt.savefig(os.path.join(self.fig_dir, f"18S.cnid.vs.num.seqs.png"), dpi=600)
        plt.close()

        # In general there are a lot of host sequences in most samples. The samples with the lowest rel abund maj seq have much fewer sequences though.
        # I think it will likely be a case that we call these warning samples and try to examine them by eye. We can also check to see if there is taxonomic
        # agreement of the sequences within these samples.

        # Another check that we will do is to look at the average pairwise distance between unique distances and see if there are any that stand out. We would
        # expect there to be a relatively low pairwise distance.
        # We should output a fasta that contains all sequences detected and we'll need to keep a track of which sequences relate to which sequence names
        # Then we can do the pairwise distance caluclation, including alignment, outside of this script and read it back in.
        # There are 9000000 unique sequences. I don't think that pairwise calculations are going to be possible.
        
        # # Dict for looking up which seq represents
        # seq_name_to_rep_seq_name_dict = {}
        # sequence_to_rep_seq_name_dict = {}
        # unique_seq_list = set()
        # # Total seq list (to make sure that all sequence names are unique). We will set it once complete and check that the lengths are the same
        # total_seq_name_list = []
        # for sample, sub_dict in sample_dict.items():
        #     print(f"Sample {sample}")
        #     for seq_name, tup in sub_dict.items():
        #         sequence = tup[3]
        #         try:
        #             # Then we have already come across this sequence
        #             rep_seq_name = sequence_to_rep_seq_name_dict[sequence]
        #             seq_name_to_rep_seq_name_dict[seq_name] = rep_seq_name
        #             total_seq_name_list.append(seq_name)
        #         except KeyError:
        #             # Then this is a new seq
        #             unique_seq_list.add(sequence)
        #             sequence_to_rep_seq_name_dict[sequence] = seq_name
        #             seq_name_to_rep_seq_name_dict[seq_name] = seq_name
        #             total_seq_name_list.append(seq_name)
        # # here we have it all populated
        # pickle.dump(seq_name_to_rep_seq_name_dict, open(os.path.join(self.resource_path, "seq_name_to_rep_seq_name_dict.p"), "wb"))
        # pickle.dump(sequence_to_rep_seq_name_dict, open(os.path.join(self.resource_path, "sequence_to_rep_seq_name_dict.p"), "wb"))
        # pickle.dump(unique_seq_list, open(os.path.join(self.resource_path, "unique_seq_list.p"), "wb"))
        # pickle.dump(total_seq_name_list, open(os.path.join(self.resource_path, "total_seq_name_list.p"), "wb"))

        # TODO a really cool figure to do would be to look at the total number of unique squences cummulatively going in order from
        # of islands and then reefs within islands to see how it increases and whether we get close to a saturation
        # We can do this by going in the correct sample order and then going through each sequene and checking to see if we've already come across it
        meta_df = pd.read_csv("/home/humebc/projects/tara/cdiv/inputs/tn_map.csv")
        meta_df = meta_df[meta_df["primers"].str.contains("18S")]
        meta_df.index = [f"TARA_{_}" for _ in list(meta_df["barcode"])]
        meta_df = meta_df.sort_values(["island", "site"], ascending = (True, True))
        # This code works but takes some time to compute
        # unique_seqs = set()
        # unique_count = 0
        # cumulative_tot = []
        # for sample in meta_df.index:
        #     print(f"{sample}")
        #     for seq, tup in sample_dict[sample].items():
        #         if tup[3] in unique_seqs:
        #             cumulative_tot.append(unique_count)
        #         else:
        #             unique_seqs.add(tup[3])
        #             unique_count += 1
        #             cumulative_tot.append(unique_count)

        # fig, ax = plt.subplots(nrows=1, ncols=1)
        # ax.plot(range(len(cumulative_tot)), cumulative_tot, 'b-', linewidth=0.5)
        # ax.set_ylabel("cummulative unique 18S seqs")
        # ax.set_xlabel("number of sequences")
        # plt.savefig(os.path.join(self.fig_dir, f"18S.novel.seqs.svg"))
        # plt.savefig(os.path.join(self.fig_dir, f"18S.novel.seqs.png"), dpi=600)
        # plt.close()


        unique_seqs = set()
        unique_count = 0
        cumulative_tot = []
        for sample in meta_df.index:
            seq = sample_to_maj_seq_dict[sample]
            if seq in unique_seqs:
                cumulative_tot.append(unique_count)
            else:
                unique_seqs.add(seq)
                unique_count += 1
                cumulative_tot.append(unique_count)

        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(range(len(cumulative_tot)), cumulative_tot, 'b-', linewidth=0.5)
        ax.set_ylabel("cummulative unique 18S maj seqs")
        ax.set_xlabel("number of sequences")
        plt.savefig(os.path.join(self.fig_dir, f"18S.novel.maj.seqs.svg"))
        plt.savefig(os.path.join(self.fig_dir, f"18S.novel.maj.seqs.png"), dpi=600)
        plt.close()

        foo = "bar"
Questions()