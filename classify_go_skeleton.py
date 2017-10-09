#!/usr/bin/python

import os
import urllib2
import argparse
import itertools


def retrieve_go_terms(protein_ids, go_terms_file="go_terms.txt"):
    """
    Get the GO data, from EBI's webservice or from the local file if it was already retrieved and saved.
    :param protein_ids: list of protein IDs.
    :param go_terms_file: file mapping uniprot ids and GO terms.
    :return: GO data read from local file or from the database.
    """

    if not go_terms_file or not os.path.isfile(go_terms_file):
        base_go_url = "http://www.ebi.ac.uk/QuickGO-Old/GAnnotation"
        go_retrieve_url = base_go_url + "?format=tsv&limit=-1&col=proteinID,goID&protein=%s"

        f = open("go_terms.txt", "w")

        for protein in protein_ids:
            response = urllib2.urlopen(go_retrieve_url % protein)
            go_data = response.read()
            print("Retrieved GO data for protein", protein)
            if go_data:
                f.write(go_data)
            else:
                print("No data found for protein", protein)

        f.close()

    else:
        with open(go_terms_file) as f:
            go_data = f.read()

    return go_data


def get_go_terms_dict(go_data):
    """
    Returns a dictionary, in which each item has a UniProt ID as its key, and
    the corresponding set of GO terms as its value. The code:
    s = go_dict["Q9BS40"]
    would store the set of GO terms associated to protein Q9BS40 in 's'.

    :param go_data: data retrieved from GO database.
    :return: dictionary with UniProt ID (key), and the corresponding set of GO terms (value).
    """
    # Define the output dictionary.
    go_dict = {}

    # Loop over the lines in 'go_data' and add every GO annotation to the appropriate set.
    for line in go_data.split('\n'):
        if not line.startswith('ID\t'):
            terms = line.split("\t")

            if len(terms) == 2:
                protein_id = terms[0]
                go_annotation = terms[1]

                # If the protein has not yet been added to the dictionary, initiate the set for the protein first.
                if protein_id not in go_dict:
                    go_dict[protein_id] = set()
                go_dict[protein_id].add(go_annotation)

    return go_dict


def read_protein_ids_file(filename):
    """
    Returns the list of UniProt IDs contained in the input file.
    :param filename: input file with UniProt IDs.
    :return: list of UniProt IDs in the input file.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    
    #######################
    ### END CODING HERE ###
    #######################
    return protein_ids


def compute_similarity_score(prot1_go_terms, prot2_go_terms):
    """
    Computes the score for two proteins on the basis of their GO terms
    :param prot1_go_terms: list of GO terms for prot1
    :param prot2_go_terms: list of GO terms for prot2
    :return: similarity score (float)
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    
    ########################
    ### END CODING HERE ####
    ########################


def check_similarity_for_protein_pair(score, threshold1, threshold2):
    """
    Given a score and two thresholds (threshold1 < threshold2), the following
    function returns the string "different" if the score is less than threshold1;
    it returns the string "ambiguous" if the score is greater than threshold1, but
    less than threshold2; it returns the string "similar" if the score is greater
    than threshold2.
    :param prot1_go_terms: list of GO terms for prot1.
    :param prot2_go_terms: list of GO terms for prot2.
    :param threshold1: bottom bound, defining different-ambiguous.
    :param threshold2: upper bound defining ambiguous-similar.
    :return: "different", "similar" or "ambiguous"

    """
    ##########################
    ### START CODING HERE ####
    ##########################

    ########################
    ### END CODING HERE ####
    ########################


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
    # You can add a pair of proteins to the list using the following code:
    # pairs_list.append((protein1, protein2))
    # Generate all possible combinations of IDs
  
    ########################
    ### END CODING HERE ####
    ########################
    return pairs


def assign_homology(go_dict, pairs, threshold1, threshold2):
    """
    Computes the similarity score between all protein pairs from the list, and decides if two proteins are homologs
    (different, ambiguous or similar).
    :param go_dict: dictionary containing the mapping of Uniprot ids (keys) and their corresponding GO terms (values)
    :param pairs: list of all possible unique protein pairs.
    :param threshold1: bottom bound, defining different-ambiguous.
    :param threshold2: upper bound defining ambiguous-similar.
    :return: dictionary with UniProt ID (key), similarity(different, ambiguous or similar) and score (value).
    """
    go_homology = {}
    ##########################
    ### START CODING HERE ####
    ##########################
  	
    ########################
    ### END CODING HERE ####
    ########################
    return go_homology


def write_results(filename, go_homology):
    """
    Writes the pairs, similarity and score to the output file.
    :param filename: name of the output file
    :param go_homology: dictionary with UniProt ID (key), similarity(different, ambiguous or similar) and score (value).
    """
    with open(filename, "w") as f:
        for pair, value in go_homology.iteritems():

            # Try/except block will process the cases when there is no data available in DB for the protein.
            # Without this block the script will stop and throw an exception for these cases.
            # For more information check: https://docs.python.org/3/tutorial/errors.html
            try:
                f.write("\t".join(pair) + "\t" + "\t".join(value) + "\n")
            except KeyError as e:
                print("No GO terms available for protein", e.args[0])


def main(input_file, output_file, go_terms_file, threshold1, threshold2):
    # Parse the input file and retrieve the GO terms associated with each protein.
    protein_ids = read_protein_ids_file(input_file)

    if go_terms_file:
        go_data = retrieve_go_terms(protein_ids, go_terms_file)
    else:
        go_data = retrieve_go_terms(protein_ids)

    go_dict = get_go_terms_dict(go_data)

    # Compute and prints the scores and the similarity string of each protein pair.
    pairs = generate_all_possible_protein_pairs(protein_ids)
    go_homology = assign_homology(go_dict, pairs, threshold1, threshold2)
    write_results(output_file, go_homology)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The script retrieves data from GO database (online or from cached file)'
                                                 ' and provides an output file'
                                                 ' with the strings "ProteinID   ProteinID   similarity  score", where'
                                                 ' similarity is a string with one of the values from '
                                                 ' different/ambiguous/similar and score is a float value.')

    parser.add_argument("-o", "--output_file", help="Output file name", required=True)
    parser.add_argument("-ids", "--protein_ids_file", help="File name of the Uniprot ids", required=True)
    parser.add_argument("-go", "--go_terms_file", help="GO database filename if it was cached", required=False)
    parser.add_argument("-t1", "--threshold1", help="Bottom bound, defining different-ambiguous", required=True)
    parser.add_argument("-t2", "--threshold2", help="Upper bound defining ambiguous-similar", required=True)

    args = parser.parse_args()

    input_file = args.protein_ids_file
    output_file = args.output_file
    go_terms_file = args.go_terms_file
    threshold1 = args.threshold1
    threshold2 = args.threshold2

    main(input_file, output_file, go_terms_file, threshold1, threshold2)