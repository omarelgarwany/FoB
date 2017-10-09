import matplotlib
matplotlib.use('AGG')

import matplotlib.pyplot as plt
import argparse
import numpy as np
import urllib2
import xml.etree.ElementTree as ET
def process_scores_file(scores_file):
	"""  
    Processing the scores file
	"""
	return [float(score) for score in open("../results/pairs_scores.csv", "r").read()[:-1].split(" ")]
def build_freq_array(scores_arr):
	"""
    Building frequency array for scores and frequencies
	"""
	freq_arr = [0 for i in range(int(max(scores_arr)*10) + 1)]
	for score in scores_arr:
		idx = int(score*10)
		freq_arr[idx] += 1
	return freq_arr
def plot(freq_arr, png_file):
	"""
    Plotting scores and frequencies
	"""
	x = np.array([float(i)/10 for i in range(11)])
	y = np.array(freq_arr)

	total_pairs = sum(y)

	y = [float(num)/total_pairs for num in y]
	print(x, y)
	plt.bar(x, y, width=0.05)
	plt.xticks(x)

	plt.xlabel("Score")
	plt.ylabel("Fraction")
    


	plt.savefig(png_file)
	return x, y
def write_to_csv(x, y, csv_file):
	"""
	Saves (score, frequency) in a CSV file
	"""
	csv_file = open(csv_file, "w")
	csv_output = ""
	for i in range(len(x)):
	    csv_output += str(x[i]) + "," + str(y[i]) + "\n"
	csv_file.write(csv_output)
	csv_file.close()

def main(scores_file, png_file, csv_file):
	scores_arr = process_scores_file(scores_file)
	freq_arr = build_freq_array(scores_arr)
	x, y = plot(freq_arr, png_file)
	write_to_csv(x, y, csv_file)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Plots frequency for each score of protein pairs')
	parser.add_argument("-scores", "--scores_file", help="Scores file path", required=True)
	parser.add_argument("-o", "--png_file", help="PNG Output file", required=True)
	parser.add_argument("-csv", "--csv_file", help="CSV File To Save Plotted Data", required=True)

	args = parser.parse_args()

	scores_file = args.scores_file
	png_file = args.png_file
	csv_file = args.csv_file
	main(scores_file, png_file, csv_file)
