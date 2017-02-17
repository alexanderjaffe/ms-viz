from pyteomics import mzxml
import json
import argparse
import glob

def retrieve_spectra(ms_run):
	
	sample = mzxml.read(ms_run)
	scans = []
	while True:
		try:
			# get next scan
			next_scan = sample.next()
			# convert np arrays to json-friendly lists
			next_scan["intensity array"] = next_scan["intensity array"].tolist()
			next_scan["m/z array"] = next_scan["m/z array"].tolist()
			scans.append(next_scan)
		except:
			break

	return(scans)

def main():

	__author__ = "Alexander L. Jaffe"
	parser = argparse.ArgumentParser(description='Converts a set of mzxml to json format for d3 viz.')
	parser.add_argument('-i','--input', help='Path to directory of mzxml.',required=True)
	parser.add_argument('-o','--output', help='Path to directory of json output.',required=False)
	args = parser.parse_args()

	if not args.output:
		out_dir = "./data"
	else:
		out_dir = args.output

	in_dir = args.input

	final_json = {}
	# select all mzXML in in_dir
	for sample in glob.glob((in_dir+"/*mzXML")):
		# get file name
		basename = sample.split("/")[(len(sample.split("/"))-1)].strip(".mzXML")
		print "Processing " + basename + ".mzXML"
		#out_path = out_dir + "/" + basename + ".json"
		# convert to json and combine into one giant file
		final_json[basename] = retrieve_spectra(sample)

	out_path = out_dir + "/" + "mzXML.json"
	with open(out_path, "w") as out_file:
		out_file.write(json.dumps(final_json))

if __name__ == '__main__':
	main()