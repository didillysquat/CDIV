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
                for ind in tax_df.index:
                    tax_string = tax_df.at[ind, "full_tax_str"]
                    try:
                        if "p_Cnidaria" in tax_string:
                            tax_dict[ind] = "Cnidaria"
                        elif "f_Symbiodiniaceae" in tax_string:
                            tax_dict[ind] = "Symbiodiniaceae"
                        else:
                            tax_dict[ind] = "other"
                    except TypeError:
                        tax_dict[ind] = "other"
                with open(os.path.join(self.post_pcr_path, f"{sample}.fasta"), "r") as f:
                    fasta_list = [_.rstrip() for _ in f]
                seq_dict = {}
                for i, line in enumerate(fasta_list):
                    if i%2 == 0:
                        seq_dict[fasta_list[i].split("\t")[0][1:]] = fasta_list[i+1]
                # Here we have a dict for each of the sequence
                for sequence in count_df.index:
                    sub_dict[sequence] = (count_dict[sequence], count_dict[sequence]/count_tot, tax_dict[sequence], seq_dict[sequence])
                sample_dict[sample] = sub_dict
            foo = "bar"
            pickle.dump(sample_dict, open(os.path.join(self.resource_path, "sample_dict.p"), "wb"))
            
        # Now we can work on the question
        # First question, how many unique sequences make up the most abudant sequence
        print("\n\n\nCounting maj seqs")
        maj_seq_dd = defaultdict(int)
        for sample, sub_dict in sample_dict.items():
            print(f"Getting maj seq for sample {sample}")
            maj_seq = [_[1][3] for _ in list(sorted(sub_dict.items(), key=lambda x: x[1][0], reverse=True)) if sub_dict[_[0]][2] == "Cnidaria"][0]
            maj_seq_dd[maj_seq] += 1

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
        foo = "bar"
        

Questions()