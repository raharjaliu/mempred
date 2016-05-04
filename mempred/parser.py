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

            if 'U' not in ss_pred:
                seq_dict[name] = (sequence, ss_pred)

        counter += 1

    return seq_dict


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
	pdbtm_sequence_dict = parse_sequence_file(pdbtm_path)

	# print opm_seq_dict
	# print pdbtm_sequence_dict

	ctr_opm = 0
	ctr_pdbtm = 0

	ctr_length_opm = 0
	ctr_length_pdbtm = 0

	opm_len_dict = {}
	pdbtm_len_dict = {}

	# variables bla and foo can be ignored. they are unimportant because they are only used as intermediate.
	bla = ""
	foo = ""

	# species is powidl. species opm dict / species pdbtm dict are the dicts for species and number of occurence.
	species = ""
	species_opm_dict = {}
	species_pdbtm_dict = {}
	for key in opm_seq_dict:
		#print key
		bla = key.split("_")[1]
		species = bla.split()[0]
		#print species
		if(species in species_opm_dict):
			species_opm_dict[species] += 1
		else:
			species_opm_dict[species] = 1
		#print species_opm_dict
		bla = opm_seq_dict[key][0]
		#print bla
		opm_len_dict[ctr_opm] = len(bla.split(",")[0])
		#bla = opm_seq_dict.values()[ctr_opm][0]
		#opm_len_dict[ctr_opm] = len(bla.split(",")[0])
		ctr_opm += 1

	#print "-----------++"

	for key2 in pdbtm_sequence_dict:
		bla = key2.split("_")[1]
		species = bla.split()[0]
		#print species
		if(species in species_pdbtm_dict):
			species_pdbtm_dict[species] += 1
		else:
			species_pdbtm_dict[species] = 1
		bla = pdbtm_sequence_dict[key2][0]
		#print bla
		# bla = pdbtm_sequence_dict.values()[ctr_pdbtm][0]
		pdbtm_len_dict[ctr_pdbtm] = len(bla.split(",")[0])
		ctr_pdbtm += 1

	"""
	Descriptive statistics for the data
	"""

	print "mean of values in opm database is " + str(np.mean(opm_len_dict.values()))
	print "median of values in opm database is " + str(np.median(opm_len_dict.values()))
	print "standard deviation of values in opm database is " + str(np.std(opm_len_dict.values()))
	print "mean of values in pdbtm database is " + str(np.mean(pdbtm_len_dict.values()))
	print "median of values in pdbtm database is " + str(np.median(pdbtm_len_dict.values()))
	print "standard deviation of values in pdbtm database is " + str(np.std(pdbtm_len_dict.values()))

	# plot histogram
	plt.hist(opm_len_dict.values(),bins=50)
	plt.title("Sequence lengths in OPM file")
	plt.xlabel("Length")
	plt.ylabel("Number of occurences")
	plt.show()

	# print species
	for i in species_pdbtm_dict:
		print str(i) + "\t" + str(species_opm_dict[i])

