#!/usr/bin/python
import os
import matplotlib.pyplot as plt
import numpy as np


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


if __name__ == "__main__":

    # path of this file's directory
    current_path = os.path.dirname(os.path.realpath(__file__))

    # path where the file opm_unmasked_hval0.fasta is located within the project
    opm_path = current_path + "/../data/sets/opm_unmasked_hval0.fasta"
    opm_buff = open(opm_path)

    # path where the file pdbtm_unmasked_hval0.fastais located within the project
    pdbtm_path = current_path + "/../data/sets/pdbtm_unmasked_hval0.fasta"
    pdbtm_buff = open(pdbtm_path)

    # parse OPM file and represent it as a map
    opm_seq_dict = parse_sequence_file(opm_path)

    # parse PDBTM file and represent it as a map
    pdbtm_seq_dict = parse_sequence_file(pdbtm_path)

    print opm_seq_dict
    print pdbtm_seq_dict


    ## Compute basic statistics for both dataset
    initial_species_statistics(opm_seq_dict, "OPM")
    initial_species_statistics(pdbtm_seq_dict, "PDBTM")


    ## Compute distribution of amino acid

    compute_sequence_distrib(opm_seq_dict, "OPM")
    compute_sequence_distrib(pdbtm_seq_dict, "PDBTM")

    print get_k_gram(opm_seq_dict[opm_seq_dict.keys()[0]][0], 3)

    ## Compute k-gram statistics of dataset

    compute_k_gram_statistics(opm_seq_dict, "OPM", 3)
    compute_k_gram_statistics(pdbtm_seq_dict, "PDBTM", 3)


    ## TODO Profile and SVM test run

    print "done"
