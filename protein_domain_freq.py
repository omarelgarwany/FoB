import matplotlib
matplotlib.use('AGG')

import matplotlib.pyplot as plt
import argparse
import numpy as np
import urllib2
import xml.etree.ElementTree as ET

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
def build_domain_count_dict(uniprot_id_filename):
    """
    Builds a dictionary containing {uniprot_id: count of different domains}
    """
    domain_count_dict = {}
    for uniprot_id in read_protein_ids_file(uniprot_id_filename):
        pfam_response_xml = retrieve_pfam_data(uniprot_id)
        domain_count_dict[uniprot_id] = len(set(parse_pfam_xml(pfam_response_xml)))
    return domain_count_dict
def build_freq_arr(uniprot_id_filename):
    """  
    Builds a frequency array for the frequency of number of unique domains 
    Example: [1, 5, 4] 
    1 protein contains 0 domains, 5 contain 1 domain, 4 contain 2 domains
    """
    domain_count_dict = build_domain_count_dict(uniprot_id_filename)
    freq_arr = [0 for i in range(max(domain_count_dict.values())+1)]
    for uniprot_id, domain_count in domain_count_dict.iteritems():
        freq_arr[domain_count] += 1
    return freq_arr
def write_to_csv(x, y, csv_file):
    """
    Writes (number of domains, frequencies) to a csv file
    """
    csv_file = open(csv_file, "w")
    csv_output = ""
    for i in range(len(x)):
        csv_output += str(x[i]) + " " + str(y[i]) + "\n"
    csv_file.write(csv_output)
    csv_file.close()
def main(uniprot_ids_file, png_file, csv_file):
    freq_arr = build_freq_arr(uniprot_ids_file)
    x = [i for i in range(len(freq_arr))]
    y = freq_arr

    plt.bar(x, y)
    plt.xticks(x)

    plt.xlabel("Number of different domains")
    plt.ylabel("Frequency")
    plt.savefig(png_file)

    write_to_csv(x, y, csv_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots frequency for each domain count of proteins')
    parser.add_argument("-ids", "--uniprot_ids_file", help="UNIPROT ids file path", required=True)
    parser.add_argument("-o", "--png_file", help="PNG Output file", required=True)
    parser.add_argument("-csv", "--csv_file", help="CSV File To Save Plotted Data", required=True)

    args = parser.parse_args()

    uniprot_ids_file = args.uniprot_ids_file
    png_file = args.png_file
    csv_file = args.csv_file
    main(uniprot_ids_file, png_file, csv_file)