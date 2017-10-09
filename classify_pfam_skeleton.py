#!/usr/bin/python

import urllib2
import xml.etree.ElementTree as ET
import argparse
import itertools


def retrieve_pfam_data(uniprot_id):
    """
    For each protein, retrieve Pfam data.
    :param uniprot_id: UniProt ID of the protein.
    :return: Pfam data corresponding to that protein in xml format.
    """

    # Constructing url to grab Pfam xml data. 
    # The format is "http://pfam.xfam.org/protein/{uniprot_id }?output=xml"
    url = "http://pfam.xfam.org/protein/" + uniprot_id.strip() + "?output=xml"

    fh = urllib2.urlopen(url)
    result = fh.read()
    fh.close()
    return result


def save_pfam_response(pfam_data, output_file):
    """
    Saves the pfam response to xml file.
    :param pfam_data: the response from the server.
    :param output_file: the output file for saving the response.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    output_file.write(pfam_data)
    return pfam_data
    ##########################
    ###  END CODING HERE  ####
    ##########################


def parse_pfam_xml(pfam_response_xml):
    """
    Parse the xml data from Pfam response.
    :param pfam_response_xml: an xml file corresponding to a protein ID extracted from Pfam database.
    :return: a dictionary containing the protein UniProtID (key) and as value a list
             containing the Pfam family/motif/domain's accession number it belongs to
             and the level of curation (Pfam-A or Pfam-B).
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    # Extracting the accession and type information from Pfam xml data 
    tree = ET.fromstring(pfam_response_xml)
    uniprot_id = tree.find("{http://pfam.xfam.org/}entry").attrib["accession"]
    matches = []
    for match in tree.iter("{http://pfam.xfam.org/}match"):
    	matches.append((match.attrib["accession"], match.attrib["type"]))
    return matches

    ##########################
    ###  END CODING HERE  ####
    ##########################





def compute_similarity_score(prot1_pfam, prot2_pfam, threshold):
    """
    Computes the score for two proteins on the basis of the data from Pfam database.
    :param prot1_pfam: data for protein 1 from Pfam.
    :param prot2_pfam: data for protein 2 from Pfam.
    :return: similarity score (float)
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    # You need to decide whether you need this function for Pfam database.
    # determines the similarity between two proteins based on the number of domains the share
    # Returns None (ambiguous) if none of the proteins have identifiabe domains (their union set is empty)
    prot1_families = set([family[0] for family in prot1_pfam])
    prot2_families = set([family[0] for family in prot2_pfam])

    

    if len(prot1_families.union(prot2_families)) == 0:
        return None
    score = float(len(prot1_families.intersection(prot2_families)))/len(prot1_families.union(prot2_families))

    return score

    ##########################
    ###  END CODING HERE  ####
    ##########################


def check_similarity_for_protein_pair(prot1_pfam, prot2_pfam, threshold):
    """
    Returns the similarity score between two proteins.
    :param pror1_pfam: pfam family/motif/domain's accession number the protein 1 belongs to
                    and the level of curation (Pfam-A or Pfam-B).
    :param pror2_pfam: pfam family/motif/domain's accession number the protein 1 belongs to
                    and the level of curation (Pfam-A or Pfam-B).
    :return: "different", "similar" or "ambiguous".
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    similarity_score = compute_similarity_score(prot1_pfam, prot2_pfam, threshold)
    # If similarity score is above threshold they are assigned a "similar" label, otherwise "different"
    # If it returns None then they are ambiguous

    if similarity_score == None:
        return "ambiguous", similarity_score
    if  similarity_score <= threshold:
        return "different", similarity_score
    elif similarity_score > threshold:
        return "similar", similarity_score
    ##########################
    ###  END CODING HERE  ####
    ##########################

    # If you will use the numeric score for Pfam (similar to GO), you may want to use check_similarity_for_protein_pair
    # with other arguments. See the example below.
    # def check_similarity_for_protein_pair(score, threshold):
    #    pass

def assign_homology(pfam_data, protein_pairs, threshold, csv_file = None):
    """
    :param pfam_data: a dictionary containing protein Uniprot IDs (keys) and the Pfam data for them(values)
    :param proteins_pairs: a file containing all protein pairs to compare
    :return: a dictionary with all protein pairs (keys) and their similarity - different/similar/ambiguous (values)
    """
    pfam_homology = dict()
    ##########################
    ### START CODING HERE ####
    ##########################
    csv_file = open(csv_file, "w")
    csv_output = ""

    for protein_pair in protein_pairs:
        similarity, similarity_score = check_similarity_for_protein_pair(pfam_data[protein_pair[0]], pfam_data[protein_pair[1]], threshold)
        ####Saving scores in a csv file to be later used to plot the scores frequencies####
        if csv_file is not None:
            if similarity_score is not None:
                csv_output += str(round(similarity_score, 1)) + " "
        ###################################################################################
        pfam_homology[protein_pair] = similarity
    csv_file.write(csv_output)
    ##########################
    ###  END CODING HERE  ####
    ##########################
    return pfam_homology


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
    for protein1_id in protein_ids:
        for protein2_id in protein_ids:
            if protein1_id != protein2_id:
                pairs.append((protein1_id, protein2_id))
    ########################
    ### END CODING HERE ####
    ########################
    return pairs


def read_protein_ids_file(filename):
    """
    Returns the list of UniProt IDs contained in the input file.
    :param filename:
    :return: list of UniProt IDs in the input file.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    protein_ids = []
    for protein_id in open(filename, "r"):
        if protein_id.endswith("\n"):
            protein_id = protein_id[0:-1]
        if protein_id.endswith(" "):
            protein_id = protein_id[0:-1] 
        protein_ids.append(protein_id)   
    #######################
    ### END CODING HERE ###
    #######################
    return protein_ids

def write_output(pfam_homology, output_file):
    """
    Writes in an output file the all of the protein pairs and their similarity/dissimilarity.
    :param pfam_homology: a dictionary with all protein pairs (keys) and
    their similarity - different/similar/ambiguous (values).
    :param output_file: the name of the output file.
    """
    with open(output_file, "w") as f:
        for pair, homology in pfam_homology.iteritems():
            f.write("\t".join([pair[0], pair[1], homology]) + "\n")




def main(protein_ids_file, output_file, xml_folder, csv_file, threshold = 0.5):
    ##########################
    ### START CODING HERE ####
    ##########################
    
    uniprot_ids = read_protein_ids_file(protein_ids_file)
    protein_lookup_dict = {}
    for uniprot_id in uniprot_ids:
        xml_data = retrieve_pfam_data(uniprot_id)
        xml_output_file = open(xml_folder + uniprot_id + ".xml", "r+")
        save_pfam_response(xml_data, xml_output_file)
        pfam_data = parse_pfam_xml(xml_data)
        protein_lookup_dict[uniprot_id] = pfam_data
    protein_pairs = generate_all_possible_protein_pairs(uniprot_ids)
    pfam_homology = assign_homology(protein_lookup_dict, protein_pairs, threshold=threshold, csv_file=csv_file)
    write_output(pfam_homology, output_file)
    #######################
    ### END CODING HERE ###
    #######################
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='The script retrieves data from Pfam database (online or from cached file)'
                                                 ' and provides an output file'
                                                 ' with the strings "ProteinID   ProteinID   similarity", where'
                                                 ' similarity is a string with one of the values from '
                                                 ' different/ambiguous/similar.')

    parser.add_argument("-ids", "--protein_ids_file", help="File name of the Uniprot ids", required=True)
    parser.add_argument("-o", "--output_file", help="Output file name", required=True)
    parser.add_argument("-xml", "--xml_folder", help="Folder where XML response from Pfam is stored", required=True)
    # We added a new argument to be able to control the similarity threshold
    parser.add_argument("-threshold", "--similarity_threshold", help="Threshold of similarity", required=False, default=0.5, type=int)
    # This argument saves scores to a csv file
    parser.add_argument("-csv", "--csv_file", help="CSV File Path", required=False, default=None)

    ##########################
    ### START CODING HERE ####
    ##########################
    # If you would like to cache the results of pfam you would probably need to add a new parameter where the cache
    # file (files) is located.
    #######################
    ### END CODING HERE ###
    #######################

    args = parser.parse_args()

    protein_ids_file = args.protein_ids_file
    output_file = args.output_file
    xml_folder = args.xml_folder
    
    threshold = args.similarity_threshold
    csv_file = args.csv_file

    main(protein_ids_file=protein_ids_file, output_file=output_file, xml_folder=xml_folder, csv_file=csv_file, threshold=threshold)

    

