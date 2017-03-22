import argparse
import json
import pandas as pd
import scipy.cluster.hierarchy as hac
import numpy as np

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

	table = pd.read_csv(cmpd_table)
	# get row and col names
	r = list(table.Compound)
	c = list(table.columns)
	
	# melt to give each observation its own entry
	melt = pd.melt(table, id_vars=["Compound", "Mass"], value_vars=list(table.columns[2:len(table.columns)]))
	# transpose
	melt["col"] = melt["Compound"].apply(lambda x: r.index(x) + 1)
	# -1 bc there are two other vars ahead of 1st in vector
	melt["row"] = melt["variable"].apply(lambda x: c.index(x) - 1)
	# small edits
	melt["var"] = melt["variable"].apply(lambda x: x.split("/")[(len(x.split("/"))-1)].replace(".mzXML",""))
	melt["cmpd"] = melt["Compound"].apply(lambda x: x.replace("compound_","#"))
	melt2 = melt.drop(["variable", "Compound"], 1)
	
	# perform hierarchical clustering on both axes
	t1 = np.matrix(table.set_index(["Compound","Mass"]))
	samp_clust = hac.linkage(t1, "average")
	t1t = np.transpose(t1)
	cmpd_clust = hac.linkage(t1t, "average")
	# add one to account for 0 indexing
	hcols = {(list(hac.leaves_list(samp_clust)).index(x) + 1):(x+1) for x in list(hac.leaves_list(samp_clust))}
	hrows = {(list(hac.leaves_list(cmpd_clust)).index(x)+1):(x+1) for x in list(hac.leaves_list(cmpd_clust))}
	# add in hclust indices - account for 0 index
	melt2["hrow"] = melt2["row"].apply(lambda x: hrows[x])
	melt2["hcol"] = melt2["col"].apply(lambda x: hcols[x])
	
	# write it out
	melt2.to_csv("../data/cmpd_table_melted.tsv", sep="\t", index=False)

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
	#processed_json = reformat_search(raw_json, keepers)
	#tabularize_search(processed_json, keepers)
	cluster_and_melt(args.input)
	print "Done."

if __name__ == '__main__':
	main()