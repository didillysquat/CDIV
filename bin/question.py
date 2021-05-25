"""
Script to answer some key quesitons about the distribution of the 18S host sequences in the CDIV data
Question 1: From the 2400 samples, roughly how many unique most abundant host 18S sequences are there
"""

import os
from numpy.lib.function_base import append
import pandas as pd
import pickle
from collections import defaultdict
import matplotlib as mpl
import matplotlib.pyplot as plt
from pandas.io.parsers import read_csv, read_table
plt.rcParams['svg.fonttype'] = 'none'
import re
from numpy import nan
import requests
from bs4 import BeautifulSoup

class QuestionsBase:
    def __init__(self):
        self.resource_path = "/home/humebc/projects/tara/cdiv/script_resources"
        self.table_output_path = "/home/humebc/projects/tara/cdiv/tables"
        self.tax_table_path = "/home/humebc/projects/tara/cdiv/output_backup_20210510/taxonomy_tables"
        self.post_pcr_path = "/home/humebc/projects/tara/cdiv/output_backup_20210510/post_pcr"
        self.fig_dir = "/home/humebc/projects/tara/cdiv/figures"
        self.sample_dict_path = os.path.join(self.resource_path, "sample_dict.p")
        self.maj_seq_dd_path = os.path.join(self.resource_path, "maj_seq_dd.p")
        self.sample_to_maj_seq_dict_path = os.path.join(self.resource_path, "sample_to_maj_seq_dict.p")
        self.sample_to_maj_seq_name_dict_path = os.path.join(self.resource_path, "sample_to_maj_seq_name_dict.p")
        self.sequence_to_tax_string_dict_path = os.path.join(self.resource_path, "sequence_to_tax_string_dict.p")
        self.seq_name_to_seq_seq_dict_path = os.path.join(self.resource_path, "seq_name_to_seq_seq_dict.p")

class Tables(QuestionsBase):
    """
    Class dedicated to making the output tables
    NB we used the Questions class to output a fasta file of the maj 18S sequences,
    along with a meta info df.
    We read this into R and produced a taxonomic group id table.
    We will work with this table as well as the pickled out objects produced in the Questions class.
    """
    def __init__(self):
        super().__init__()
        self.tax_group_df = read_csv(os.path.join(self.resource_path, "tax.id.groups.csv"))
        self.tax_group_df = self.tax_group_df.iloc[:,1:]
        self.tax_group_df = self.tax_group_df.set_index("tax_group_ID")

        # Now we want to create the other tables
        if (
            os.path.exists(self.maj_seq_dd_path) and 
            os.path.exists(self.sample_to_maj_seq_dict_path) and 
            os.path.exists(self.sample_to_maj_seq_name_dict_path) and 
            os.path.exists(self.sequence_to_tax_string_dict_path) and 
            os.path.exists(self.sample_dict_path) and
            os.path.exists(self.seq_name_to_seq_seq_dict_path)
            ):
            # This is the key dictionary that holds the sequence abundances and their meta information for every sample.
            # The other dictionaries are derivatives from this dictionary
            # They key is a sample-id, the value is a dict. This dict is keys of sequence name(from the fastq) to a tuple of:
            # absolute seq abundance; relative seq abundance; seq tax assignation as Cnidaria, Symbiodiniaceae or other; the nucleotide sequence; the annotation string from mmseqs;
            self.sample_dict = pickle.load(open(self.sample_dict_path, "rb"))
            # This is k= nucleotide sequences, v= int number of samples the seq was maj Cnidarian seq in
            self.maj_seq_dd = pickle.load(open(self.maj_seq_dd_path, "rb"))
            # k = sample-id v=nucleotide seq of most abund cnidarian seq
            self.sample_to_maj_seq_dict = pickle.load(open(self.sample_to_maj_seq_dict_path, "rb"))
            # k = sample-id v=fastq name of most abund cnidarian seq
            self.sample_to_maj_seq_name_dict = pickle.load(open(self.sample_to_maj_seq_name_dict_path, "rb"))
            # k = nucleotide seq to 
            self.sequence_to_tax_string_dict = pickle.load(open(self.sequence_to_tax_string_dict_path, "rb"))
            # We will generate a futher dictionary that is the seq_name to the seq_seq from the self.tax_group_df
            self.seq_name_to_seq_seq_dict = pickle.load(open(self.seq_name_to_seq_seq_dict_path, "rb"))
            # We also want the inverse of this dict
            self.seq_seq_to_seq_name_dict = {v: k for k, v in self.seq_name_to_seq_seq_dict.items()}
            # We also want to be able to pick out which tax id group a given sequence is in
            # This will allow us to associate the above info to the self.tax_group_df
            self.seq_name_to_tax_group_dict = {}
            for tax_group, ser in self.tax_group_df.iterrows():
                list_of_seq_names = ser["maj_18S_seq_list_names"].split(",")
                for seq_name in list_of_seq_names:
                    self.seq_name_to_tax_group_dict[seq_name] = tax_group
            # TODO get the picture URLs from somewhere.
            
        else:
            raise RuntimeError("The pickled resources do not exist. You will need to run the Questions class to create them.")

        # To be able to look up the photos we need to know the sampling-design_label
        # We can look this up in the provenance table.
        # Once we have the sampling-design label, we can then look up the sample_label for the pictures that are associated with the coral sample
        # We can then build the URL to the pictures from that.
        # It would be good if we could check that the URL at least exists.
        # Pangea site that holds the photos is https://store.pangaea.de/Projects/TARA-PACIFIC/Images/
        # We will beautiful soup that link and get all of the .jpg files available to download from that page and use these as to verify the existnce of the URLs.
        jpg_links_pickle_path = os.path.join(self.resource_path, "jpg_links.p")
        self.pangea_pic_base_URL = "https://store.pangaea.de/Projects/TARA-PACIFIC/Images/"
        if os.path.exists(jpg_links_pickle_path):
            jpg_links = pickle.load(open(jpg_links_pickle_path, "rb"))
        else:
            soup = BeautifulSoup(requests.get(self.pangea_pic_base_URL).text, features="html.parser")
            links = [_.string for _ in soup.find_all('a') if _.string]
            jpg_links = [str(_) for _ in links if _.endswith(".jpg")]
            pickle.dump(jpg_links, open(jpg_links_pickle_path, "wb"))

        # Load in the provenance table
        provenance_df_slimmed_pickle_path = os.path.join(self.resource_path, "provenance_df_slimmed.p")
        if os.path.exists(provenance_df_slimmed_pickle_path):
            self.provenance_df_slimmed = pickle.load(open(provenance_df_slimmed_pickle_path, "rb"))
        else:
            self.provenance_df = read_table("/home/humebc/projects/tara/cdiv/inputs/TARA-PACIFIC_samples-provenance_20200731d.txt", skiprows=[0])
            # We will slim down this provenance df to just the sample-ids and associated image files that we are working with and then pickle out
            # To do this we will need to get a list of the sampling-design labels that correspond to the sample-id and then slim according to these.
            # This is because the associated pictures will have the same sample-design label as the samples they correspond to.
            sample_list = list(self.sample_to_maj_seq_dict.keys())
            sampling_design_labels = list(self.provenance_df[self.provenance_df["sample-id_source"].isin(sample_list)]["sampling-design_label"].values)
            self.provenance_df_slimmed = self.provenance_df[self.provenance_df["sampling-design_label"].isin(sampling_design_labels)]
            self.provenance_df_slimmed.set_index("sample-id_source", drop=False, inplace=True)
            pickle.dump(self.provenance_df_slimmed, open(provenance_df_slimmed_pickle_path, "wb"))       
        
        # Now we can make a sample-id to URL list dict
        sample_to_pic_url_list_dict_pickle_path = os.path.join(self.resource_path, "sample_to_pic_url_list_dist.p")
        if os.path.exists(sample_to_pic_url_list_dict_pickle_path):
            self.sample_to_pic_url_list_dict = pickle.load(open(sample_to_pic_url_list_dict_pickle_path, "rb"))
        else:
            self.sample_to_pic_url_list_dict = defaultdict(list)
            missing_pictures = []
            samples_with_no_pictures = []
            salvaged_samples = []
            for sample in self.sample_to_maj_seq_dict.keys():
                sampling_design_label = self.provenance_df_slimmed.at[sample, "sampling-design_label"]
                non_sample_rows = self.provenance_df_slimmed[(self.provenance_df_slimmed["sampling-design_label"] == sampling_design_label) & (self.provenance_df_slimmed["sample-id_source"] != sample)]
                # It may be that there are no pictures in which case the above would be empty
                if len(non_sample_rows.index) == 0:
                    # In this case we will manually search through the jpg list to see if we find matching pictures because it could be that pictures have been uploaded but not
                    # assigned sample names in the provenance table or something similar.
                    # NB this salvaged 13 samples.
                    jpegs_with_name = [_ for _ in jpg_links if sampling_design_label in _]
                    if jpegs_with_name:
                        salvaged_samples.append(sample)
                        self.sample_to_pic_url_list_dict[sample] = [os.path.join(self.pangea_pic_base_URL, _) for _ in jpegs_with_name]
                    else:
                        samples_with_no_pictures.append(sample)
                        
                # The picture names will be the sample_label plus .jpg
                for pic_sample, ser in non_sample_rows.iterrows():
                    # Check to see if the sample pic exists
                    pic_name = ser["sample_label"] + ".jpg"
                    if pic_name in jpg_links:
                        self.sample_to_pic_url_list_dict[sample].append(os.path.join(self.pangea_pic_base_URL, pic_name))
                    else:
                        missing_pictures.append(os.path.join(self.pangea_pic_base_URL, pic_name))

        foo = "bar"

        # Add the picture URLs to the self.tax_group_df
        # TODO. Also when we do the self.tax_group_df URLs, check to see how many of the tax groups don't have picture URLs. If we get lucky, it may be 0.
        # Go tax group by tax group.
        # Get the list of samples and generate the picture URL list from this.
        # also generate data for the following columns: number_samples_with_pictures	samples_with_pictures	number_samples_without_pictures	samples_without_pictures	picture_URLs
        number_samples_with_pictures_list = []
        samples_with_pictures_list = []
        number_samples_without_pictures_list = []
        samples_without_pictures_list = []
        picture_URLs_list = []
        for tax_group, ser in self.tax_group_df.iterrows():
            sample_list = ser["sample_list"].split(",")
            samples_without_pictures = [_ for _ in sample_list if _ not in self.sample_to_pic_url_list_dict.keys()]
            samples_without_pictures_list.append(",".join(samples_without_pictures))
            number_samples_without_pictures_list.append(len(samples_without_pictures))
            samples_with_pictures = [_ for _ in sample_list if _ in self.sample_to_pic_url_list_dict.keys()]
            samples_with_pictures_list.append(",".join(samples_with_pictures))
            number_samples_with_pictures_list.append(len(samples_with_pictures))
            picture_URLs_list.append(",".join([",".join(self.sample_to_pic_url_list_dict[_]) for _ in samples_with_pictures]))
        
        self.tax_group_df["number_samples_with_pictures"] = number_samples_with_pictures_list
        self.tax_group_df["samples_with_pictures"] = samples_with_pictures_list
        self.tax_group_df["number_samples_without_pictures"] = number_samples_without_pictures_list
        self.tax_group_df["samples_without_pictures"] = samples_without_pictures_list
        self.tax_group_df["picture_URLs"] = picture_URLs_list
        # The tax_group_df is now complete so output
        self.tax_group_df.to_csv(os.path.join(self.table_output_path, "taxonomy_groups_18S.tsv"), sep="\t")

        # We want to produce two futher tables
        # Table one is the 'by sample' table
        # Table two is the 'by sequence' table
        
        # By sample table
        
        sample_id_list = []
        maj_cnid_18S_seq_seq_list = []
        maj_cnid_18S_seq_name_list = []
        mmseq_tax_annot_string_list = []        
        bioinf_tax_annot_phylum_list = []
        bioinf_tax_annot_class_list = []
        bioinf_tax_annot_order_list = []
        bioinf_tax_annot_family_list = []
        bioinf_tax_annot_genus_list = []
        bioinf_tax_annot_species_list = []
        bioinf_lowest_rank_list = []
        tax_group_ID_list = []
        evidence_for_tax_ID_group_list = []
        has_been_visually_IDed_list = []
        visual_ID_based_on_sample_id_list = []
        visual_tax_phylum_list = []    
        visual_tax_class_list = []
        visual_tax_order_list = []
        visual_tax_family_list = []
        visual_tax_genus_list = []
        visual_tax_species_list = []
        visual_lowest_tax_rank_list = []
        has_pictures_list = []
        picture_URLs_list = []            
        for sample_id, seq in self.sample_to_maj_seq_dict.items():
            sample_id_list.append(sample_id)
            maj_cnid_18S_seq_seq_list.append(seq)
            maj_cnid_18S_name = self.seq_seq_to_seq_name_dict[seq]
            maj_cnid_18S_seq_name_list.append(maj_cnid_18S_name)
            mmseq_tax_annot_string = self.sequence_to_tax_string_dict[seq]
            mmseq_tax_annot_string_list.append(mmseq_tax_annot_string)
            # Here we need to extract the phylumn, class, order, family, genus and species if they are available
            # otherwise we ill put an na.
            (
                bioinf_tax_annot_phylum, bioinf_tax_annot_class, bioinf_tax_annot_order,
                bioinf_tax_annot_family, bioinf_tax_annot_genus, bioinf_tax_annot_species
                ) = [self._pull_out_tax_from_tax_string(reg_ex, mmseq_tax_annot_string) for reg_ex in ["(p_[A-Za-z]+);", "(c_[A-Za-z]+);", "(o_[A-Za-z]+);", "(f_[A-Za-z]+);", "(g_[A-Za-z]+);", "(s_[A-Za-z]+ [A-Za-z]+)"]]
            bioinf_tax_annot_phylum_list.append(bioinf_tax_annot_phylum)
            bioinf_tax_annot_class_list.append(bioinf_tax_annot_class)
            bioinf_tax_annot_order_list.append(bioinf_tax_annot_order)
            bioinf_tax_annot_family_list.append(bioinf_tax_annot_family)
            bioinf_tax_annot_genus_list.append(bioinf_tax_annot_genus)
            bioinf_tax_annot_species_list.append(bioinf_tax_annot_species)
            
            tax_annotations = [bioinf_tax_annot_phylum, bioinf_tax_annot_class, bioinf_tax_annot_order, bioinf_tax_annot_family, bioinf_tax_annot_genus, bioinf_tax_annot_species]
            tax_level_dict = {0: "phylum", 1: "class", 2: "order", 3: "family", 4: "genus", 5:"species"}
            for i, annot in enumerate(tax_annotations):
                lowest_tax = 0
                if pd.isna(annot):
                    # Then the lowest level was the previous tax level
                    bioinf_lowest_rank = tax_level_dict[i-1]
                    break
                bioinf_lowest_rank = "species"
            bioinf_lowest_rank_list.append(bioinf_lowest_rank)
            
            tax_group_ID = self.seq_name_to_tax_group_dict[maj_cnid_18S_name]
            tax_group_ID_list.append(tax_group_ID)
            evidence_for_tax_ID_group_list.append(self.tax_group_df.at[tax_group_ID, "evidence_for_tax_ID_group"])

            has_been_visually_IDed_list.append("no")
            visual_ID_based_on_sample_id_list.append(nan)
            (visual_tax_phylum, visual_tax_class, visual_tax_order, visual_tax_family, visual_tax_genus, visual_tax_species, visual_lowest_tax_rank) = [nan for _ in range(7)]
            visual_tax_phylum_list.append(visual_tax_phylum)
            visual_tax_class_list.append(visual_tax_class)
            visual_tax_order_list.append(visual_tax_order)
            visual_tax_family_list.append(visual_tax_family)
            visual_tax_genus_list.append(visual_tax_genus)
            visual_tax_species_list.append(visual_tax_species)
            visual_lowest_tax_rank_list.append(visual_lowest_tax_rank)
            
            has_pictures = False
            picture_URLs = nan
            if sample_id in self.sample_to_pic_url_list_dict.keys():
                has_pictures = True
                picture_URLs = ",".join(self.sample_to_pic_url_list_dict[sample_id])
            has_pictures_list.append(has_pictures)
            picture_URLs_list.append(picture_URLs)

        self.by_sample_df = pd.DataFrame()
        self.by_sample_df ["sample-id"] = sample_id_list
        self.by_sample_df["maj_cnid_18S_seq_seq"] = maj_cnid_18S_seq_seq_list
        self.by_sample_df["maj_cnid_18S_seq_name"] = maj_cnid_18S_seq_name_list
        self.by_sample_df["mmseq_tax_annot_string"] = mmseq_tax_annot_string_list
        self.by_sample_df["bioinf_tax_annot_phylum"] = bioinf_tax_annot_phylum_list
        self.by_sample_df["bioinf_tax_annot_class"] = bioinf_tax_annot_class_list
        self.by_sample_df["bioinf_tax_annot_order"] = bioinf_tax_annot_order_list
        self.by_sample_df["bioinf_tax_annot_family"] = bioinf_tax_annot_family_list
        self.by_sample_df["bioinf_tax_annot_genus"] = bioinf_tax_annot_genus_list
        self.by_sample_df["bioinf_tax_annot_species"] = bioinf_tax_annot_species_list
        self.by_sample_df["bioinf_lowest_rank"] = bioinf_lowest_rank_list
        self.by_sample_df["tax_group_ID"] = tax_group_ID_list
        self.by_sample_df["evidence_for_tax_ID_group"] = evidence_for_tax_ID_group_list
        self.by_sample_df["has_been_visually_IDed"] = has_been_visually_IDed_list
        self.by_sample_df["visual_ID_based_on_sample-id"] = visual_ID_based_on_sample_id_list
        self.by_sample_df["visual_tax_phylum"] = visual_tax_phylum_list
        self.by_sample_df["visual_tax_class"] = visual_tax_class_list
        self.by_sample_df["visual_tax_order"] = visual_tax_order_list
        self.by_sample_df["visual_tax_family"] = visual_tax_family_list
        self.by_sample_df["visual_tax_genus"] = visual_tax_genus_list
        self.by_sample_df["visual_tax_species"] = visual_tax_species_list
        self.by_sample_df["visual_tax_lowest_rank"] = visual_lowest_tax_rank_list
        self.by_sample_df["has_pictures"] = has_pictures_list
        self.by_sample_df["picture_URLs"] = picture_URLs_list
        self.by_sample_df.set_index("sample-id", inplace=True, drop=True)

        # output the by sample df
        self.by_sample_df.to_csv(os.path.join(self.table_output_path, "by_sample_18S.tsv"), sep="\t", index=True)
        foo = "bar"

        # TODO maybe get rid of the sequence name as this is just a bit confusing. We can always look up the sequence name for our own purposes.

    def _pull_out_tax_from_tax_string(self, reg_ex_pattern, mmseq_tax_annot_string):
        match_obj = re.search(reg_ex_pattern, mmseq_tax_annot_string)
        if match_obj:
            group_obj = match_obj.groups()
            if len(group_obj) > 1:
                raise RuntimeError("group contains more than 1 match")
            else:
                # Remove the "p_" part and return
                return group_obj[0][2:]
        else:
            return nan


class Questions(QuestionsBase):
    def __init__(self):
        super().__init__()

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

        # We will also want to make a phylogenetic tree based on the maj 18S sequences
        # We will output a fasta that has the abundance of each of the samples built into its name
        # We will also put out a meta info df to work with in R where we will visualize the tree
        # and generate the taxonomic group tables.
        # We want to add a list of samples that each sequence represents
        # We will make a dict to do that
        seq_seq_to_sample_list_dict = defaultdict(list)
        seq_name_to_seq_seq_dict = {}
        for k, v in sample_to_maj_seq_dict.items():
            seq_seq_to_sample_list_dict[v].append(k)
        
        fasta_list = []
        seq_count = 0
        tree_df_data = [[],[],[],[],[]]
        for seq, abund in sorted_maj_seqs:
            seq_name_to_seq_seq_dict[f"seq{seq_count}_{abund}"] = seq
            fasta_list.append(f">seq{seq_count}_{abund}")
            fasta_list.append(f"{seq}")
            tree_df_data[0].append(f"seq{seq_count}_{abund}")
            tree_df_data[1].append(seq)
            tree_df_data[2].append(abund)
            tree_df_data[3].append(sequence_to_tax_string_dict[seq])
            tree_df_data[4].append(",".join(seq_seq_to_sample_list_dict[seq]))
            seq_count += 1
        # Pickle out the seq_name_to_seq_seq_dict
        pickle.dump(seq_name_to_seq_seq_dict, open(self.seq_name_to_seq_seq_dict_path, "wb"))
        with open(os.path.join(self.resource_path, "18S.maj.seqs.fasta"), "w") as f:
            for line in fasta_list:
                f.write(f"{line}\n")
        meta_df = pd.DataFrame(list(zip(*tree_df_data[1:])), columns=["seq", "abund", "tax_string", "sample_list"], index=tree_df_data[0])
        meta_df.to_csv(os.path.join(self.resource_path, "18S.maj.seqs.meta.df.csv"))
        
        # We'd like to know what the average proportion the most abundant sequence represents in a sample
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

Tables()