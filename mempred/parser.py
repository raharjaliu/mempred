import os


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

    print opm_seq_dict
    print pdbtm_sequence_dict


    """
    Both OPM and PDBTM files are now parsed and represented as dictionary (very similar to Map in Java).
    The format of the dictionary is SEQUENCE_NAME -> (SEQUENCE, SECONDARY_STRUCTURE_PREDICTION)

    TODO -- Quirin

    DO FOLLOWING STATS:

    - number of sequence in each of the databases
    - length statistics: MEAN, MEDIAN, STANDARD DEVIATION
    - histogram of length distribution and save it:
        read http://matplotlib.org/1.2.1/examples/pylab_examples/histogram_demo.html to learn how to plot histogram in python
    - how many occurrences of each species:
        PDB species can be seen by the name. For example, following sequence
            >G1UBD1|3rgbI|PMOB_METCA 0.97
        belongs to the species METCA. Use function some_string.split() to split a given list

    If done, upload the statistics and histogram plot onto group wiki page: https://i12r-studfilesrv.informatik.tu-muenchen.de/pp12016bio/index.php/Group_2
    """
