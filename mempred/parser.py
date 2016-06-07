#!/usr/bin/python
import os
import matplotlib.pyplot as plt
import numpy as np
import random as rd


def parse_sequence_file(file_path):
    """
    Given a secondary structure FASTA file, parses the file
    and returns the contents of it as dictionary.
    The dictionary is defined as mapping from
    sequence header to tupple of sequence and its secondary
    structure prediction.
    """

    in_buff = open(file_path)
    in_content = in_buff.readlines()

    counter = 0
    seq_dict = {}

    for line in in_content:

        # remove tailing new line
        line = line.strip()

        if counter % 3 == 0:

            name = line

        elif counter % 3 == 1:

            sequence = line

        elif counter % 3 == 2:

            ss_pred = line

            seq_dict[name] = (sequence, ss_pred)

            """
            if 'U' not in ss_pred:
                seq_dict[name] = (sequence, ss_pred)
            """

        counter += 1

    return seq_dict


def get_k_gram(strin, k):
    """
    Compute k-gram of a given sequence
    """

    if k > len(strin):

        return {}

    else:

        kgram = {}

        for i in xrange(len(strin) - k + 1):

            this_str = strin[i:(i + k)]

            if this_str not in kgram:

                kgram[this_str] = 1

            else:

                kgram[this_str] += 1

        return kgram



def get_aa_distrib(seq):
    """
    Calculate the distribution of
    amino acid in a given sequence
    """

    aa_distrib = {'a': 0, 'r': 0, 'n': 0,
                  'd': 0, 'c': 0, 'e': 0,
                  'q': 0, 'g': 0, 'h': 0,
                  'i': 0, 'l': 0, 'k': 0,
                  'm': 0, 'f': 0, 'p': 0,
                  's': 0, 't': 0, 'w': 0,
                  'y': 0, 'v': 0}

    for chr in seq:

        aa_distrib[chr] = aa_distrib[chr] + 1

    return aa_distrib

def merge_non_uniform_dict(d1, d2):

    """
    Merge the values of two dicts.
    Unlike merge_dict, it is not assumed
    that both dicts have identical key set
    """

    keys = set(d1.keys() + d2.keys())

    merged_dict = {}

    for key in keys:

        if key in d1:

            merged_dict[key] = d1[key]

        if key in d2:

            if key in merged_dict:

                merged_dict[key] = merged_dict[key] + d2[key]

            else:

                merged_dict[key] = d2[key]

    return merged_dict


def merge_dict(d1, d2):

    """
    Merge the values of two dicts.
    It is strictly assumed that both dicts have identical keyset
    and the values are of numerical types
    """

    for key in d2:

        d1[key] = d1[key] + d2[key]

    return d1

def initial_species_statistics(seq_dict, dataset_name):
    """
    Get initial statistics for a dataset.

    Data set is represented as dictionary with following structure:

    NAME : (SEQUENCE, SECONDARY_STRUCTURE_ANNOTATION)
    """

    print "=" * 10
    print "STATS FOR %s" % dataset_name

    counter = 0
    len_dict = {}
    species_dict = {}

    for key in seq_dict:
        # print key
        species = (key.split("_")[1]).split()[0]
        # print species
        if species in species_dict:
            species_dict[species] += 1
        else:
            species_dict[species] = 1

        this_seq = seq_dict[key][0]
        len_dict[counter] = len(this_seq)
        counter += 1

    print "MEAN: " + str(np.mean(len_dict.values()))
    print "MEDIAN: " + str(np.median(len_dict.values()))
    print "STD: " + str(np.std(len_dict.values()))

    plt.hist(len_dict.values(), bins=50)
    plt.title("Sequence Length Histogram")
    plt.xlabel("Fequency")
    plt.ylabel("Length")
    # plt.show()
    plt.savefig("length_hist_%s.png" % dataset_name, bbox_inches='tight')
    plt.close()

    for key in species_dict:
        print "%s\t%s" % (key, str(species_dict[key]))

    print "=" * 10


def compute_sequence_distrib(seq_dict, dataset_name):
    """
    Compute basic sequence statistics
    for a given dataset. The statistics computed are:
    - mean sequence length
    - median sequence length
    - std of sequence length

    The function also computes the frequency of
    UniProt species to appear in the dataset.
    """

    seq_distrib = None

    for seq_name in seq_dict.keys():

        this_seq = seq_dict[seq_name][0]
        this_distrib = get_aa_distrib(this_seq.lower())

        if seq_distrib is None:

            seq_distrib = this_distrib

        else:

            seq_distrib = merge_dict(seq_distrib, this_distrib)

    keys = []
    vals = []

    for key in seq_distrib.keys():
        keys += [key]
        vals += [seq_distrib[key]]

    plt.bar(range(len(vals)), vals)
    plt.xticks(range(len(vals)), keys)
    plt.title("Amino Acid Distribution of %s Dataset" % dataset_name)
    plt.xlabel("Amino Acid")
    plt.ylabel("Frequency")
    plt.savefig("amino_acid_distrib_%s.png" % dataset_name, bbox_inches='tight')
    plt.close()


def compute_k_gram_statistics(seq_dict, dataset_name, k):
    """
    Compute k-gram statistics of data set
    """

    seq_k_gram = None

    for seq_name in seq_dict:

        this_seq = seq_dict[seq_name][0]
        this_k_gram = get_k_gram(this_seq, 3)

        if seq_k_gram is None:

            seq_k_gram = this_k_gram

        else:

            seq_k_gram = merge_non_uniform_dict(seq_k_gram, this_k_gram)

    keys = []
    vals = []

    for key in seq_k_gram.keys():
        keys += [key]
        vals += [seq_k_gram[key]]

    plt.hist(vals, bins=len(vals))

    # plt.bar(range(len(vals)), vals)
    # plt.xticks(range(len(vals)), [])
    plt.title("Histogram of Sequence %d-Gram for Dataset %s" % (k, dataset_name))
    plt.ylabel("Frequency")
    # plt.xlabel("k-Gram")
    plt.savefig("%d_gram_of_%s.png" % (k, dataset_name), bbox_inches='tight')
    plt.close()


def create_initial_profile_ranges(ss_annot):
    """
    Given secondary structure annotation sequence,
    compute initial profile ranges. Initial profile
    ranges describe the starting and end point of
    Helical ("positive") and non-Helical ("negaitve")
    annotations

    E.g. for secondary structure annotation of

    HHHHH11111UUUUU11111HHHHH

    the positive ranges would be [(0, 4), (20, 24)]
    while the negtive ranges would be [(5, 19)]
    """

    pos_index = []
    neg_index = []
    is_pos = False
    last_index = 0

    for index in range(len(ss_annot)):

        if index == 0:
            is_pos = ss_annot[index] == 'H'

        if ss_annot[index] == 'H':

            if not is_pos:

                is_pos = True

                if index > 0:
                    neg_index += [(last_index, index - 1)]
                    last_index = index

        else:

            if is_pos:

                is_pos = False

                if index > 0:
                    pos_index += [(last_index, index - 1)]
                    last_index = index

    if is_pos:

        pos_index += [(last_index, len(ss_annot) - 1)]

    else:

        neg_index += [(last_index, len(ss_annot) - 1)]

    return pos_index, neg_index


def extract_legal_profile_ranges(init_profile_range, profile_length=20):
    """
    Given intial profile ranges from function

    create_initial_profile_range()

    create profile ranges that conform with
    project requirement
    """

    extracted_index = []

    for this_start, this_end in init_profile_range:

        length = this_end - this_start + 1

        if length >= profile_length:

            this_index = this_start

            sample_num = length / profile_length
            # the total flexibility of sampling
            flexibility = length % profile_length
            # the single fliexibility range for each sampling
            single_flexibility = flexibility / sample_num

            for i in range(sample_num):
                this_margin = rd.randint(0, single_flexibility)
                this_sample_start = this_index + this_margin
                this_sample_end = this_sample_start + profile_length - 1

                extracted_index += [(this_sample_start, this_sample_end)]

                this_index = this_sample_end + 1


    return extracted_index

def compute_preliminiary_statistics(dataset, dataset_name, k_gram=3):

    ## Compute basic statistics for both dataset
    initial_species_statistics(dataset, dataset_name)

    ## Compute distribution of amino acid
    compute_sequence_distrib(dataset, dataset_name)

    ## Compute k-gram statistics of dataset
    compute_k_gram_statistics(dataset, dataset_name, k_gram)


def write_profiles(profiles_path, output_path, dataset_name, sequence_dict, part_training=.9):

    file_list = os.listdir(profiles_path)
    file_list = [profiles_path + this_file for this_file in file_list if 'blastPsiMat' in this_file]

    pos_profiles = {}
    neg_profiles = {}

    pos_seq = {}
    neg_seq = {}

    for this_file in file_list:

        this_file_seq_id = this_file.split('/')[-1].split('.')[0].split('_')[1]
        (this_seq, this_ss) = sequence_dict[this_file_seq_id]

        (this_pos_index, this_neg_index) = create_initial_profile_ranges(this_ss)


        this_pos_index = extract_legal_profile_ranges(this_pos_index, profile_length=20)
        this_neg_index = extract_legal_profile_ranges(this_neg_index, profile_length=20)

        # print "this_ss: %s" % this_ss
        # print this_pos_index
        # print this_neg_index
        # print "=" * 10

        in_buff = open(this_file)
        in_content = in_buff.readlines()

        header = "".join(in_content[0:3])
        footer = "".join(in_content[-6:])

        content_lines = in_content[3:(-6)]

        for this_index in this_pos_index:

            # print "#" * 10
            # print this_index
            # print "#" * 10

            # print this_index
            this_content = "".join(content_lines[this_index[0]:(this_index[1] + 1)])

            this_profile = header + this_content + footer
            # this_profile = header + this_content
            # this_seq = sequence_dict[this_file_seq_id][0][(this_index[0] - 1):this_index[1]]
            this_seq = sequence_dict[this_file_seq_id][0][this_index[0]:(this_index[1] + 1)]

            # print this_profile
            # print this_seq

            pos_profiles[this_file_seq_id + "_" + str(this_index[0]) + "_" + str(this_index[1])] = this_profile
            pos_seq[this_file_seq_id + "_" + str(this_index[0]) + "_" + str(this_index[1])] = this_seq

        for this_index in this_neg_index:

            # print "#" * 10
            # print this_index
            # print "#" * 10

            # print this_index
            this_content = "".join(content_lines[this_index[0]:(this_index[1] + 1)])

            this_profile = header + this_content + footer
            # this_profile = header + this_content
            # this_seq = sequence_dict[this_file_seq_id][0][(this_index[0] - 1):this_index[1]]
            this_seq = sequence_dict[this_file_seq_id][0][this_index[0]:(this_index[1] + 1)]

            # print this_profile
            # print this_seq

            neg_profiles[this_file_seq_id + "_" + str(this_index[0]) + "_" + str(this_index[1])] = this_profile
            neg_seq[this_file_seq_id + "_" + str(this_index[0]) + "_" + str(this_index[1])] = this_seq

        # print '=' * 20

    len_dev_set_pos = int(part_training * len(pos_profiles.keys()))
    len_dev_set_neg = int(part_training * len(neg_profiles.keys()))

    ## NEW TEST

    dev_seqs = pos_profiles.keys()[:len_dev_set_pos] + neg_profiles.keys()[:len_dev_set_neg]
    rd.shuffle(dev_seqs)
    test_seqs = pos_profiles.keys()[len_dev_set_pos:] + neg_profiles.keys()[len_dev_set_neg:]
    rd.shuffle(test_seqs)

    fileout_string_dev = ""
    fileout_string_dev_list = ""
    arff_dev_id = "@attribute identifier {"
    arff_dev_dat = "@data"
    fileout_string_dev_arff = ""

    fileout_string_test = ""
    fileout_string_test_list = ""
    arff_test_id = "@attribute identifier {"
    arff_test_dat = "@data"
    fileout_string_test_arff = ""

    for seq_name in dev_seqs:

        arff_dev_id += (seq_name + ", ")
        fileout_string_dev_list += (">" + seq_name + "\n")

        if seq_name in pos_profiles.keys():

            fileout_string_dev += (">" + seq_name + "\n" + pos_seq[seq_name] + "\n" + pos_profiles[seq_name])
            arff_dev_dat += ("\n" + seq_name + ",positive")

        else:

            fileout_string_dev += (">" + seq_name + "\n" + neg_seq[seq_name] + "\n" + neg_profiles[seq_name])
            arff_dev_dat += ("\n" + seq_name + ",negative")

    fileout_string_dev_list = fileout_string_dev_list[:-1]
    arff_dev_id = arff_dev_id[:-2] + "}"
    fileout_string_dev_arff = "@relation docs\n" + arff_dev_id + "\n@attribute class {positive,negative}\n" + arff_dev_dat

    text_file = open(output_path + dataset_name + "_dev.profile", "w")
    text_file.write(fileout_string_dev)
    text_file.close()

    text_file = open(output_path + dataset_name + "_dev.arff", "w")
    text_file.write(fileout_string_dev_arff)
    text_file.close()

    text_file = open(output_path + dataset_name + "_dev.list", "w")
    text_file.write(fileout_string_dev_list)
    text_file.close()

    for seq_name in test_seqs:

        arff_test_id += (seq_name + ", ")
        fileout_string_test_list += (">" + seq_name + "\n")

        if seq_name in pos_profiles.keys():

            fileout_string_test += (">" + seq_name + "\n" + pos_seq[seq_name] + "\n" + pos_profiles[seq_name])
            arff_test_dat += ("\n" + seq_name + ",positive")

        else:

            fileout_string_test += (">" + seq_name + "\n" + neg_seq[seq_name] + "\n" + neg_profiles[seq_name])
            arff_test_dat += ("\n" + seq_name + ",negative")

    fileout_string_test_list = fileout_string_test_list[:-1]
    arff_test_id = arff_test_id[:-2] + "}"
    fileout_string_test_arff = "@relation docs\n" + arff_test_id + "\n@attribute class {positive,negative}\n" + arff_test_dat

    text_file = open(output_path + dataset_name + "_test.profile", "w")
    text_file.write(fileout_string_test)
    text_file.close()

    text_file = open(output_path + dataset_name + "_test.arff", "w")
    text_file.write(fileout_string_test_arff)
    text_file.close()

    text_file = open(output_path + dataset_name + "_test.list", "w")
    text_file.write(fileout_string_test_list)
    text_file.close()

    print "DATASET: %s" % dataset_name
    print "NUMBER OF SEQ IN POSITIVE TRAINING SET: %d" % len_dev_set_pos
    print "NUMBER OF SEQ IN NEGATIVE TRAINING SET: %d" % len_dev_set_neg
    print "NUMBER OF SEQ IN POSITIVE TESTING SET: %d" % (len(pos_profiles.keys()) - len_dev_set_pos)
    print "NUMBER OF SEQ IN NEGATIVE TESTING SET: %d" % (len(neg_profiles.keys()) - len_dev_set_neg)
    print "TOTAL NUMBER OF SEQ IN TRAINING SET: %d" % (len_dev_set_pos + len_dev_set_neg)
    print "TOTAL NUMBER OF SEQ IN TESTING SET: %d" % (len(pos_profiles.keys()) + len(neg_profiles.keys()) - len_dev_set_pos - len_dev_set_neg)


    ## END OF NEW TEST

    # fileout_string_dev_pos = ""
    # fileout_string_dev_pos_list = ""
    # fileout_string_test_pos = ""
    # fileout_string_test_pos_list = ""
    # fileout_string_dev_neg = ""
    # fileout_string_dev_neg_list = ""
    # fileout_string_test_neg = ""
    # fileout_string_test_neg_list = ""
    #
    # for seq_name in pos_profiles.keys()[:len_dev_set_pos]:
    #     this_string = ">" + seq_name + "\n" + pos_seq[seq_name] + "\n" + pos_profiles[seq_name]
    #     fileout_string_dev_pos += this_string
    #     fileout_string_dev_pos_list += (">" + seq_name + "\n")
    #
    # for seq_name in pos_profiles.keys()[len_dev_set_pos:]:
    #     this_string = ">" + seq_name + "\n" + pos_seq[seq_name] + "\n" + pos_profiles[seq_name]
    #     fileout_string_test_pos += this_string
    #     fileout_string_test_pos_list += (">" + seq_name + "\n")
    #
    # for seq_name in neg_profiles.keys()[:len_dev_set_neg]:
    #     this_string = ">" + seq_name + "\n" + neg_seq[seq_name] + "\n" + neg_profiles[seq_name]
    #     fileout_string_dev_neg += this_string
    #     fileout_string_dev_neg_list += (">" + seq_name + "\n")
    #
    # for seq_name in neg_profiles.keys()[len_dev_set_neg:]:
    #     this_string = ">" + seq_name + "\n" + neg_seq[seq_name] + "\n" + neg_profiles[seq_name]
    #     fileout_string_test_neg += this_string
    #     fileout_string_test_neg_list += (">" + seq_name + "\n")
    #
    #
    # print "DATASET: %s" % dataset_name
    # print "NUMBER OF SEQ IN POSITIVE TRAINING SET: %d" % len_dev_set_pos
    # print "NUMBER OF SEQ IN POSITIVE TESTING SET: %d" % (len(pos_profiles.keys()) - len_dev_set_pos)
    # print "NUMBER OF SEQ IN NEGATIVE TRAINING SET: %d" % len_dev_set_neg
    # print "NUMBER OF SEQ IN NEGATIVE TESTING SET: %d" % (len(neg_profiles.keys()) - len_dev_set_neg)
    #
    # text_file = open(output_path + dataset_name + "_pos_dev.profile", "w")
    # text_file.write(fileout_string_dev_pos[:-1])
    # text_file.close()
    #
    # text_file = open(output_path + dataset_name + "_pos_dev.list", "w")
    # text_file.write(fileout_string_dev_pos_list[:-1])
    # text_file.close()
    #
    # text_file = open(output_path + dataset_name + "_pos_test.profile", "w")
    # text_file.write(fileout_string_dev_pos[:-1])
    # text_file.close()
    #
    # text_file = open(output_path + dataset_name + "_pos_test.list", "w")
    # text_file.write(fileout_string_dev_pos_list[:-1])
    # text_file.close()
    #
    # text_file = open(output_path + dataset_name + "_neg_dev.profile", "w")
    # text_file.write(fileout_string_dev_neg[:-1])
    # text_file.close()
    #
    # text_file = open(output_path + dataset_name + "_neg_dev.list", "w")
    # text_file.write(fileout_string_dev_neg_list[:-1])
    # text_file.close()
    #
    # text_file = open(output_path + dataset_name + "_neg_test.profile", "w")
    # text_file.write(fileout_string_dev_neg[:-1])
    # text_file.close()
    #
    # text_file = open(output_path + dataset_name + "_neg_test.list", "w")
    # text_file.write(fileout_string_dev_neg_list[:-1])
    # text_file.close()


if __name__ == "__main__":

    # path of this file's directory
    current_path = os.path.dirname(os.path.realpath(__file__))

    # path where the file opm_unmasked_hval0.fasta is located within the project
    opm_path = current_path + "/../data/sets/opm_unmasked_hval0.fasta"

    # path where the file pdbtm_unmasked_hval0.fastais located within the project
    pdbtm_path = current_path + "/../data/sets/pdbtm_unmasked_hval0.fasta"

    # parse OPM file and represent it as a map
    opm_seq_dict = parse_sequence_file(opm_path)

    # parse PDBTM file and represent it as a map
    pdbtm_seq_dict = parse_sequence_file(pdbtm_path)

    # print opm_seq_dict
    # print pdbtm_seq_dict

    ## Compute basic statistics for both dataset
    # compute_preliminiary_statistics(opm_seq_dict, "OPM", k_gram=3)
    # compute_preliminiary_statistics(pdbtm_seq_dict, "PDBTM", k_gram=3)


    ## Divide both dataset into training and test sets

    # delete existing file in dataset folder

    dataset_path = current_path + "/../dataset/"
    file_list = os.listdir(dataset_path)
    file_list = [dataset_path + this_file for this_file in file_list]

    for this_file in file_list:

        if os.path.isfile(this_file):
            os.unlink(this_file)


    # create modified dictionary of sequences in which the key refers true sequence name

    opm_dict_real_key = {}

    for key in opm_seq_dict:

        new_key = key.split('|')[0][1:]
        opm_dict_real_key[new_key] = opm_seq_dict[key]

    pdbtm_dict_real_key = {}

    for key in pdbtm_seq_dict:

        new_key = key.split('|')[0][1:]
        pdbtm_dict_real_key[new_key] = pdbtm_seq_dict[key]

    # print opm_dict_real_key
    # print pdbtm_dict_real_key

    ## create profile dict from original profile maps

    part_training_pdbtm = .9
    profiles_path = current_path + "/../TM_proteins_profiles/"
    output_path = dataset_path
    dataset_name = "PDBTM"

    write_profiles(profiles_path, dataset_path, "PDBTM", pdbtm_dict_real_key, part_training=0.9)
    write_profiles(profiles_path, dataset_path, "OPM", opm_dict_real_key, part_training=0.9)

    print "done"
