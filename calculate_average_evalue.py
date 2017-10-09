#This script calculates the average of an evalue distribution
import numpy as np
import argparse

def calculate_average(csv_file):
	print("Calculating")
	return np.mean([float(a) for a in open(csv_file, "r").readlines()])
def main(csv_file):
    print("Average of evalues is "  + str(calculate_average(csv_file)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates average for a group of evalues")
    parser.add_argument("-csv", "--csv_file", help="CSV File that contains all evalues", required=True)

    args = parser.parse_args()
    csv_file = args.csv_file
    main(csv_file)
