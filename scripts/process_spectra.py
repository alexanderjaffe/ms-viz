'''PROCESSES SPECTRA FOR VISUALIZATION/GNPS'''

from pyteomics import mzxml
import json
import argparse
import glob
import sys
import math
import os

# given cmpd table and spectra file,
# generates keyed list of spectra to keep
def get_keepers(cmpd_table, spectra_file):

   	# read in cmpd table
   	with open(cmpd_table) as f:
   		table = f.readlines()
   	f.close()

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
					#sample = path.split("/")[len(path.split("/"))-1].replace(".mzXML",'')
					sample = os.path.basename(path)
					# take first spectrum as representative
					spectrum = element.split("$")[1].split(",")[0]
					# build sample: spectra to keep dict
					if sample not in keepers.keys():
						keepers[sample] = {spectrum: cmpd}
					else:
						keepers[sample][spectrum] = cmpd
	f2.close()

	return(keepers)

# given keepers and sample names,
# actually retrieves spectra to keep
def retrieve_spectra(ms_run, basename, keepers_list):

	sample = mzxml.read(ms_run)
	scans = []

	#print keepers_list[basename].keys()

	# read through mzXML file
	while True:
		try:
			# get next scan
			next_scan = sample.next()
			# filter for keepers
			if next_scan["num"] in keepers_list[basename].keys():

				# convert np arrays to json-friendly lists, and round vals
				ia = [round(val,2) for val in next_scan["intensity array"]]
				mz = [math.floor(val*10)/10 for val in next_scan["m/z array"]]
				# write peaks as json pairs
				if len(ia) == len(mz):
					spectrum = {}
					for i in range(len(mz)):
						# only write out pair if int. > 0
						if (ia[i]>0):
							# sum intensities from close mzs into combined peaks
							if mz[i] not in spectrum:
								spectrum[mz[i]] = ia[i]
							else:
								spectrum[mz[i]] += ia[i]
					# write out as json pairs
					spectrum_list = []
					for key, item in spectrum.iteritems():
						spectrum_list.append({"i":item, "mz":key})

					next_scan["spectrum"] = spectrum_list
					next_scan["compound"] = keepers_list[basename][next_scan["num"]]
					next_scan["sample"] = basename
					next_scan.pop("intensity array")
					next_scan.pop("m/z array")

					scans.append(next_scan)
		except:
			break

	return(scans)

# convert list of processed spectra to mgf
def spectra_as_mgf(spectra, basename):

	outfile = open("../data/" + basename + ".mgf","w")
	i = 0
	for spectrum in spectra:
		if i != 0: # separate entries
			outfile.write("\n")
		# write spectrum block
		outfile.write("BEGIN IONS\n")
		pre = spectrum["precursorMz"][0]
		outfile.write("PEPMASS=%s %s\n" %(pre["precursorMz"], pre["precursorIntensity"]))
		try:
			outfile.write("CHARGE=%s%s\n" %(pre["precursorCharge"], spectrum["polarity"]))
		except:
			outfile.write("CHARGE=%s%s\n" %("?", spectrum["polarity"]))
		outfile.write("MSLEVEL=2\n")
		outfile.write("FILENAME=%s\n" %(spectrum["sample"]))
		outfile.write("ACTIVATION=%s\n" %(pre["activationMethod"]))
		outfile.write("INSTRUMENT=orbitrap\n")
		outfile.write("SCANS=%s\n" %(spectrum["num"]))
		for peak in spectrum["spectrum"]:
			outfile.write(str(peak["mz"]) + " " + str(peak["i"]) + "\n")
		outfile.write("END IONS\n")
		i += 1
	outfile.close()

def main():

	__author__ = "Alexander L. Jaffe"
	parser = argparse.ArgumentParser(description='Converts a set of mzxml to json format for d3 viz and mgf files for gnps querying.')
	parser.add_argument('-i','--input', help='Path to directory of mzxml.',required=True)
	parser.add_argument('-c','--compound_table', help='Path to filtered compound table.',required=True)
	parser.add_argument('-s','--spectra_file', help='Path to spectra mapping file.',required=True)

	args = parser.parse_args()
	in_dir = args.input

	# generate keepers list for each sample
	keepers_dict = get_keepers(args.compound_table, args.spectra_file)

	final_json = []
	# select all mzXML in in_dir
	for sample in glob.glob((in_dir+"/*mzXML")):
		# get file name
		#basename = sample.split("/")[(len(sample.split("/"))-1)].replace(".mzXML", '')
		basename = os.path.basename(sample)
		print "Processing " + basename
		# convert to json and combine into one giant file
		results = retrieve_spectra(sample, basename, keepers_dict)
		final_json = final_json + results
		print "Total of %d spectra..." %(len(final_json))

		# for each sample, create mgf for gnps analysis
		spectra_as_mgf(results,basename)

	# write out spectra in json form
	out_path = "../data/spectra.json"
	with open(out_path, "w") as out_file:
		out_file.write(json.dumps(final_json))
	out_file.close()

if __name__ == '__main__':
	main()
