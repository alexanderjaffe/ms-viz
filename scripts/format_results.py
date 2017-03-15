import argparse
import json
import pandas as pd

def reformat_search(raw_json, keepers):

	json_out = open("../data/gnps_LS_processed.json","w")

	final_json = []
	for line in raw_json:
		temp = {}
		# select key/value pairs to keep
		for keys,item in line.iteritems():
			if keys in keepers:
				temp[keys] = item
		final_json.append(temp)

	# write out processed json
	json_out.write(json.dumps(final_json))
	json_out.close()

	return final_json

def tabularize_search(processed_json, fields):

	table_out = open("../data/gnps_LS.table","w")
	
	# write out header
	headers = "\t".join(fields)
	table_out.write(headers + "\n")

	# write out lines of table
	for entry in processed_json:
		temp = []
		for field in fields:
			try:
				temp.append(entry[field])
			except:
				temp.append(None)
		line = "\t".join(temp)
		table_out.write(line + "\n")

	table_out.close()

# perform hierarchical clustering on cmpd table
# then reformat for visualization input
def cluster_and_melt(cmpd_table):

	table = pd.read


def main():

	__author__ = "Alexander L. Jaffe"
	parser = argparse.ArgumentParser(description='Reformats GNPS and compound table files for visualization.')
	parser.add_argument('-i','--input', help='Path to filtered compound table.',required=True)
	args = parser.parse_args()

	in_file = open("../data/gnps_LS_raw.json")
	raw_json = json.load(in_file)

	# define fields to keep from gnps response json
	keepers = ["id", "SpectrumFile", "#Scan#", "Compound_Name", "LibraryQualityString", \
		"MQScore", "SharedPeaks", "TIC_Query", "SpecMZ", "LibMZ", "MassDiff", "Charge"]
	
	print "Processing library search and compound table for visualization..."
	processed_json = reformat_search(raw_json, keepers)
	tabularize_search(processed_json, keepers)
	cluster_and_melt(args.input)
	print "Done."

if __name__ == '__main__':
	main()