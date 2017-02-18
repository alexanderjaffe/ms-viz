from pyteomics import mzxml
import json
import argparse
import glob
import sys

def retrieve_spectra(ms_run, basename, keepers_list):
	
	sample = mzxml.read(ms_run)
	scans = []
	while True:
		try:
			# get next scan
			next_scan = sample.next()
			# convert np arrays to json-friendly lists
			next_scan["intensity array"] = next_scan["intensity array"].tolist()
			next_scan["m/z array"] = next_scan["m/z array"].tolist()
			# filter for keepers
			if next_scan["num"] in keepers_dict[basename]:
				scans.append(next_scan)
			else: pass 
		except:
			break

	return(scans)

def get_keepers(cmpd_table, spectra_file):

	'''# help from alex CC - preprocess_data.py
	f = open(cmpd_table)
	i=0
	cmpds, firstline = [],[]
	for line in f.readlines():
		if i!=0:
			cmpds.append(line.split(",")[0])
		else:
			samples.append(line.split(","))
			i+=1

	samples = []
	for item in firstline:
		if "mzXML" in item:
			samples.append(item.split["/"][len(item.split["/"]-1)].strip(".mzXML"))
	print cmpds, samples '''

   	with open(cmpd_table) as f:
   		table = f.readlines()

   	# check for filtered cmpd_table
   	if len(table) > 200:
   		print "Error: Make sure to use compound table filtered for size/singletons."
   		sys.exit()

   	# get cmpds to keep from cmpd_table
   	cmpds = []
   	for row in table:
   		if "mzXML" not in row:
   			cmpds.append(row.split(",")[0])

	# get spectra to keep from spectra_file
	keepers = {}
	f2 = open(spectra_file)
	for line in f2.readlines():
		cmpd = line.split("|")[0]
		# only keep cmpds from cmpd table
		if cmpd in cmpds:
			for element in line.strip("\n").split("|"):
				if "mzXML" in element:
					path = element.split("$")[0]
					sample = path.split("/")[len(path.split("/"))-1].strip(".mzXML")
					spectra = element.split("$")[1].split(",")
					# build sample: spectra to keep dict
					if sample not in keepers.keys():
						keepers[sample] = spectra
					else:
						keepers[sample] = keepers[sample] + spectra

	return(keepers)

def main():

	__author__ = "Alexander L. Jaffe"
	parser = argparse.ArgumentParser(description='Converts a set of mzxml to json format for d3 viz.')
	parser.add_argument('-i','--input', help='Path to directory of mzxml.',required=True)
	parser.add_argument('-c','--compound_table', help='Path to filtered compound table.',required=True)
	parser.add_argument('-s','--spectra_file', help='Path to spectra mapping file.',required=True)
	#parser.add_argument('-o','--output', help='Path to output json.',required=False)

	args = parser.parse_args()

	'''if not args.output:
		out_dir = "./data"
	else:
		out_dir = args.output'''

	in_dir = args.input

	# generate keepers list for each sample
	keepers_dict = get_keepers(args.compound_table, args.spectra_file)

	final_json = {}
	# select all mzXML in in_dir
	for sample in glob.glob((in_dir+"/*mzXML")):
		# get file name
		basename = sample.split("/")[(len(sample.split("/"))-1)].strip(".mzXML")
		print "Processing " + basename + ".mzXML"
		#out_path = out_dir + "/" + basename + ".json"
		# convert to json and combine into one giant file
		#final_json[basename] = retrieve_spectra(sample, basename, keepers_dict)

	#out_path = out_dir + "/" + "mzXML.json"
	out_path = "mzXML.json"
	#with open(out_path, "w") as out_file:
		#out_file.write(json.dumps(final_json))

if __name__ == '__main__':
	main()