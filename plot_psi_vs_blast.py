import matplotlib
matplotlib.use('AGG')

import matplotlib.pyplot as plt
import argparse
import numpy as np

#This script produces a combined ROC plot for both blast and psi-blast given a points files
def plot_psi_vs_blast(psi_points, blast_points, png_filename):
    """
    This functions plots both BLAST and PSI-BLAST ROC plots
    """
    x_blast, y_blast = blast_points
    x_psi, y_psi = psi_points

    plt.plot(x_blast, y_blast, label="BLAST")
    plt.plot(x_psi, y_psi, label="PSI-BLAST")

    plt.plot([0,1],[0,1],'--k')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    plt.legend()
    plt.savefig(png_filename)

def process_points_file(points_file):
    """
    Processes the blast and psi blast points file (A file with FPR [x points] vs TPR [y points] for different evalues)
    """
    x, y = [a.split(" ") for a in open(points_file, "r").readlines()]
    return [float(num) for num in x], [float(num) for num in y]

def main(psi_points_file, blast_points_file, png_filename):
    psi_points = process_points_file(psi_points_file)
    blast_points = process_points_file(blast_points_file)
    plot_psi_vs_blast(psi_points, blast_points, png_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw and save a ROC plot to a file for both blast and psi-blast files together")
    parser.add_argument("-psi", "--psi_points_file", help="PSI-BLAST points file", required=True)
    parser.add_argument("-blast", "--blast_points_file", help="BLAST points file", required=True)
    parser.add_argument("-o", "--png_filename", help="PNG File Name", required=True)

    args = parser.parse_args()

    psi_points_file = args.psi_points_file
    blast_points_file = args.blast_points_file
    png_filename = args.png_filename

    main(psi_points_file, blast_points_file, png_filename)