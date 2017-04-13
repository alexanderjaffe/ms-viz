'''FORMATS GNPS RESULTS AND CMPD TABLE FOR VISUALIZATION'''

import argparse
import json
import pandas as pd
import scipy.cluster.hierarchy as hac
import numpy as np

# given raw json, and a list of fields,
# generates cleaned json
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

# given processed json and fields, writes out table
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
def cluster_and_melt(cmpd_table, min_mz, max_mz):

	intab = pd.read_csv(cmpd_table)
	
	# filtering based on mz
	# convert mass to number
	intab["temp"] = intab["Mass"].apply(lambda x: float(x.split("+-")[0]))
	temp = intab[(intab.temp > float(min_mz)) & (intab.temp < float(max_mz))]
	table = temp.drop("temp", 1)

	# get row and col names
	r = list(table.Compound)
	c = list(table.columns)

	# compute col indices for most 'conserved' cmpds
	cmpd_dicts = []
	for keys, items in table.iterrows():
		temp_dict = {"cmpd":items.Compound, "sample_count":0}
	
		# if > 0, cmpd is 'present'
		for i in range(2,(len(items)-1)):
			if items[i] > 0:
				temp_dict["sample_count"] +=1

		cmpd_dicts.append(temp_dict)

	# now sort cmpds by number of times they appear across samples
	cmpd_dicts_sorted = sorted(cmpd_dicts, key=lambda k: k['sample_count'], reverse=True) 
	cmpd_list_sorted = [item["cmpd"] for item in cmpd_dicts_sorted]
	# get 'ranking' for a compound based on how it sorted
	table["ccol"] = table["Compound"].apply(lambda x: (cmpd_list_sorted.index(x)+1))
	
	# melt to give each observation its own entry
	melt = pd.melt(table, id_vars=["Compound", "Mass","ccol"])
	# transpose
	melt["col"] = melt["Compound"].apply(lambda x: r.index(x) + 1)
	# -1 bc there are two other vars ahead of 1st in vector
	melt["row"] = melt["variable"].apply(lambda x: c.index(x) - 1)
	# small edits
	melt["var"] = melt["variable"].apply(lambda x: x.split("/")[(len(x.split("/"))-1)].replace(".mzXML",""))
	melt["cmpd"] = melt["Compound"].apply(lambda x: x.replace("compound_","#"))
	melt["log_value"] = melt["value"].apply(lambda x: np.log10((x+1)))
	melt2 = melt.drop(["variable", "Compound"], 1)
	
	# perform hierarchical clustering on both axes
	t1 = np.matrix(table.set_index(["Compound","Mass","ccol"]))
	samp_clust = hac.linkage(t1, "ward")
	t1t = np.transpose(t1)
	cmpd_clust = hac.linkage(t1t, "ward")
	# add one to account for 0 indexing
	hcols = {(x+1):(list(hac.leaves_list(samp_clust)).index(x) + 1) for x in list(hac.leaves_list(samp_clust))}
	hrows = {(x+1):(list(hac.leaves_list(cmpd_clust)).index(x) +1) for x in list(hac.leaves_list(cmpd_clust))}
	
	# add in hclust indices - account for 0 index
	melt2["hrow"] = melt2["row"].apply(lambda x: hrows[x])
	melt2["hcol"] = melt2["col"].apply(lambda x: hcols[x])
	
	# write it out
	melt2.to_csv("../data/cmpd_table_melted.tsv", sep="\t", index=False)

def main():

	__author__ = "Alexander L. Jaffe"
	parser = argparse.ArgumentParser(description='Reformats GNPS and compound table files for visualization.')
	parser.add_argument('-i','--input', help='Path to filtered compound table.',required=True)
	parser.add_argument('-min','--minimum_mz', help='Minimum mz to keep.',required=False)
	parser.add_argument('-max','--maximum_mz', help='Maximum mz to keep.',required=False)
	args = parser.parse_args()

	in_file = open("../data/gnps_LS_raw.json")
	raw_json = json.load(in_file)

	# define fields to keep from gnps response json
	keepers = ["id", "SpectrumFile", "#Scan#", "Compound_Name", "LibraryQualityString", \
		"MQScore", "SharedPeaks", "TIC_Query", "SpecMZ", "LibMZ", "MassDiff", "Charge"]
	
	# process filtering args
	# if none supplied, use really high/low values
	min_mz = args.minimum_mz if args.minimum_mz else 0
	max_mz = args.maximum_mz if args.maximum_mz else 1e10

	print "Processing library search and compound table for visualization..."
	# run the functions
	processed_json = reformat_search(raw_json, keepers)
	tabularize_search(processed_json, keepers)
	cluster_and_melt(args.input, min_mz, max_mz)
	print "Done."

if __name__ == '__main__':
	main()