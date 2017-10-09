#!/usr/bin/python

import itertools
import argparse


def retrieve_scop_data(scop_file):
    """
    Reads a databaset file from SCOP and returns a dictionary with the protein IDS mapping.
    :param scop_file: database file containing mapping of PDB's to SCOP ID's.
    :return: dictionary containing the mapping of PDBs (keys) and their corresponding SCOP data (values).
    """

    scop_data = {}

    ##########################
    ### START CODING HERE ####
    ##########################
    # You can parse SCOP data in various ways. E.g. you can use dictionary of dictionaries
    # {proteinID: {"class": class, "fold": fold, "superfamily": superfamily, 'family': family}}

    ########################
    ### END CODING HERE ####
    ########################

    return scop_data


def compute_similarity_score(prot1_scop, prot2_scop):
    """
    Computes the score for two proteins on the basis of the data from SCOP database.
    :param prot1_scop: data for protein 1 from SCOP database.
    :param prot2_scop: data for protein 2 from SCOP database.
    :return: similarity score (float)
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    # You need to decide whether you need this function for SCOP database.
    pass

    ########################
    ### END CODING HERE ####
    ########################


def check_similarity_for_protein_pair(prot1_scop, prot2_scop):
    """
    Returns the similarity score between two proteins.
    :param prot1_scop: data for protein 1 from SCOP database.
    :param prot2_scop: data for protein 2 from SCOP database.
    :param pair: a tuple with the UniProt IDs of the two proteins to compare.
    :return: "different", "similar" or "ambiguous".
    """
    ##########################
    ### START CODING HERE ####
    ##########################

    ########################
    ### END CODING HERE ####
    ########################

# If you will use the numeric score for SCOP (similar to GO), you may want to use check_similarity_for_protein_pair
# with other arguments. See the example below.
# def check_similarity_for_protein_pair(score, threshold):
#    pass

def generate_all_possible_protein_pairs(protein_ids):
    """
    Returns a list containing all unique protein pairs.
    :param protein_ids: list of all proteins IDs.
    :return: list of possible unique protein pairs.
    """
    pairs = list()
    ##########################
    ### START CODING HERE ####
    ##########################

    ########################
    ### END CODING HERE ####
    ########################
    return pairs


def assign_homology(scop_dict, protein_ids_pdbs, pairs):
    """
    Computes the similarity score between all protein pairs from the list, and decides if two proteins are homologs
    (different, ambiguous or similar).
    :param scop_dict: dictionary containing the mapping of PDBs (keys) and their corresponding SCOP data (values).
    :param protein_ids_pdbs: dictionary with UniprotID as key and PDB ID as a value.
    :param pairs: list of all possible unique protein pairs.
    :return: dictionary with UniProt ID (key), similarity(different, ambiguous or similar).
    """
    scop_homology = {}

    ##########################
    ### START CODING HERE ####
    ##########################
    # You should remember to take care about the proteins that are not in the SCOP database.
 
    ########################
    ### END CODING HERE ####
    ########################
    return scop_homology


def write_results(filename, scop_homology):
    """
    Writes in an output file the all of the protein pairs and their similarity/dissimilarity.
    :param output_file: the name of the output file.
    :param scop_homology: dictionary (keys: protein pairs as tuples; values: one of the value - different/similar/ambiguous)
    """
    with open(filename, "w") as f:
        for (p1, p2), value in scop_homology.iteritems():
            f.write("\t".join([p1, p2, value]) + "\n")


def read_protein_ids_file(filename):
    """
    Returns the list of UniProt IDs contained in the input file.
    :param filename:
    :return: list of UniProt IDs in the input file.
    """
    ##########################
    ### START CODING HERE ####
    ##########################

    #######################
    ### END CODING HERE ###
    #######################
    return protein_ids


def read_lookup_table(filename):
    """
    Reads the specified file and returns the dictionary with UniprotID as key and PDB ID as a value.
    :param filename: file with the mapping between Uniprot ids and PDB ids.
    :return: dictionary with UniprotID as key and PDB ID as a value.
    """
    ##########################
    ### START CODING HERE ####
    ##########################

    #######################
    ### END CODING HERE ###
    #######################
    return protid_pdb


def main(input_file, output_file, pdb_id_file, scop_file):
    ##########################
    ### START CODING HERE ####
    ##########################

    #######################
    ### END CODING HERE ###
    #######################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The script retrieves data from SCOP database (from local file)'
                                                 ' and provides an output file'
                                                 ' with the strings "ProteinID   ProteinID   similarity", where'
                                                 ' similarity is a string with one of the values from '
                                                 ' different/ambiguous/similar.')

    parser.add_argument("-o", "--output_file", help="Output file name")
    parser.add_argument("-ids", "--protein_ids_file", help="File with the protein Uniprot IDs")
    parser.add_argument("-pdb", "--pdb_id_file", help="File with the mapping between Uniprot ids and PDB ids")
    parser.add_argument("-s", "--scop_file", help="SCOP database file")

    args = parser.parse_args()

    input_file = args.protein_ids_file
    output_file = args.output_file
    pdb_id_file = args.pdb_id_file
    scop_file = args.scop_file

    main(input_file, output_file, pdb_id_file, scop_file)