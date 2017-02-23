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
			# filter for keepers
			if next_scan["num"] in keepers_list[basename].keys():
				# convert np arrays to json-friendly lists, and round vals
				ia = [round(val,2) for val in next_scan["intensity array"]]
				mz = [round(val,2) for val in next_scan["m/z array"]]
				# write peaks as json pairs
				if len(ia) == len(mz):
					spectrum = []
					for i in range(len(mz)):
						spectrum.append({"i":ia[i], "mz":mz[i]})
					next_scan["spectrum"] = spectrum
					next_scan["compound"] = keepers_list[basename][next_scan["num"]]
					next_scan["sample"] = basename
					next_scan.pop("intensity array")
					next_scan.pop("m/z array")
				scans.append(next_scan)
		except:
			break
	return(scans)

def get_keepers(cmpd_table, spectra_file):

   	with open(cmpd_table) as f:
   		table = f.readlines()
   	f.close()

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
			# iterate through spectra entries and add to dict 
			for element in line.strip("\n").split("|"):
				if "mzXML" in element:
					path = element.split("$")[0]
					# parse out sample name and spectra
					sample = path.split("/")[len(path.split("/"))-1].strip(".mzXML")
					# take first spectrum as representative
					spectrum = element.split("$")[1].split(",")[0]
					# build sample: spectra to keep dict
					if sample not in keepers.keys():
						keepers[sample] = {spectrum: cmpd}
					else:
						keepers[sample][spectrum] = cmpd
	f2.close()
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

	final_json = []
	# select all mzXML in in_dir
	for sample in glob.glob((in_dir+"/*mzXML")):
		# get file name
		basename = sample.split("/")[(len(sample.split("/"))-1)].strip(".mzXML")
		print "Processing " + basename + ".mzXML"
		#out_path = out_dir + "/" + basename + ".json"
		# convert to json and combine into one giant file
		#final_json[basename] = retrieve_spectra(sample, basename, keepers_dict)
		final_json = final_json + retrieve_spectra(sample, basename, keepers_dict)
		print "Total of %d spectra..." %(len(final_json))

	#out_path = out_dir + "/" + "mzXML.json"
	out_path = "../data/spectra.json"
	with open(out_path, "w") as out_file:
		out_file.write(json.dumps(final_json))
	out_file.close()
	
if __name__ == '__main__':
	main()