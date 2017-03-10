import ftplib
import argparse
import sys
import requests
import json
import time
import getpass
import os
import glob

def ftp(mode, path, basename, login, password):

	ftp = ftplib.FTP("ccms-ftp01.ucsd.edu")
	try:
		ftp.login(login, password)
	except:
		print "Login failed :/"
		sys.exit()
	if mode == "upload":
		myfile = open(path, "rb")
		try:
			ftp.storbinary("STOR " + basename, myfile)
		except:
			print "File upload failed :/"
			sys.exit()
		myfile.close()
	elif mode == "delete":
		ftp.delete(basename)
	ftp.quit()

# invokes GNPS workflow
# BY M. WANG - miw023@cs.ucsd.edu
def invoke_workflow(base_url, parameters, login, password):
	username = login
	password = password

	# turn off verification warnings
	requests.packages.urllib3.disable_warnings()
	s = requests.Session()

	payload = {
		'user' : login,
		'password' : password,
		'login' : 'Sign in'
	}

	r = s.post('https://' + base_url + '/ProteoSAFe/user/login.jsp', data=payload, verify=False)
	r = s.post('https://' + base_url + '/ProteoSAFe/InvokeTools', data=parameters, verify=False)
	task_id = r.text

	if len(task_id) > 4 and len(task_id) < 60:
		#print("Launched Task: : " + r.text)
		return task_id
	else:
		print(task_id)
		return None

# creates parameter mapping for GNPS workflow
# ADAPTED FROM M. WANG - miw023@cs.ucsd.edu
def launch_workflow(username, file_name, workflow):
	
	# set general gnps params
	parameters_map = {}

	parameters_map["ANALOG_SEARCH"] = 0
	parameters_map["email"] = "alexander_jaffe@berkeley.edu"
	parameters_map["desc"] = "Query spectra for d3 viz."
	parameters_map["FILTER_LIBRARY"] = 0
	parameters_map["FILTER_PRECURSOR_WINDOW"] = 1
	parameters_map["FILTER_STDDEV_PEAK_INT"] = 0
	parameters_map["library_on_server"] = "d.speclibs"
	parameters_map["MAX_SHIFT_MASS"] = 100.0
	parameters_map["MIN_MATCHED_PEAKS"] = 6
	parameters_map["MIN_PEAK_INT"] = 50.0
	parameters_map["SCORE_THRESHOLD"] = 0.7
	parameters_map["spec_on_server"] = "f." + username + "/" + file_name
	parameters_map["tolerance.Ion_tolerance"] = 0.5
	parameters_map["tolerance.PM_tolerance"] = 1.5
	parameters_map["uuid"] = "1DAB2BC6-6827-0001-9AB5-390F1E781419"
	parameters_map["WINDOW_FILTER"] = 1
	#parameters_map["protocol"] = "none"
	#parameters_map["TASK_ID_COMPARISON"] = task_comparison
	#parameters_map["USER_ID_COMPARISON"] = "continuous"

	# set workflow specific params
	if workflow == "search":
		parameters_map["SEARCH_LIBQUALITY"] = 3
		parameters_map["TOP_K_RESULTS"] = 1
		parameters_map["workflow"] = "MOLECULAR-LIBRARYSEARCH"
	if workflow == "network":
		parameters_map["CLUSTER_MIN_SIZE"] = 2
		parameters_map["CREATE_CLUSTER_BUCKETS"] = 0
		parameters_map["CREATE_TOPOLOGY_SIGNATURES"] = 0
		parameters_map["FIND_MATCHES_IN_PUBLIC_DATA"] = 0
		parameters_map["MAXIMUM_COMPONENT_SIZE"] = 0
		parameters_map["MIN_MATCHED_PEAKS_SEARCH"] = 3
		parameters_map["PAIRS_MIN_COSINE"] = 0.80
		parameters_map["RUN_MSCLUSTER"] = "off"
		parameters_map["TOP_K"] = 1
		parameters_map["SEARCH_LIBQUALITY"] = 3
		parameters_map["workflow"] = "METABOLOMICS-SNETS"

	return(parameters_map)

# Waits for not running, then returns status
# BY M. WANG - miw023@cs.ucsd.edu
def wait_for_workflow_finish(base_url, task_id):
	url = 'https://' + base_url + '/ProteoSAFe/status_json.jsp?task=' + task_id
	json_obj = json.loads(requests.get(url, verify=False).text)
	while (json_obj["status"] != "FAILED" and json_obj["status"] != "DONE"):
		print("Working...")
		time.sleep(10)
		try:
			json_obj = json.loads(requests.get(url, verify=False).text)
		except KeyboardInterrupt:
			raise
		except:
			print("Exception In Wait")
			time.sleep(1)

	return json_obj["status"]

def main():

	__author__ = "Alexander L. Jaffe"
	parser = argparse.ArgumentParser(description='Performs GNPS library search or molecular networking on spectra files.')
	parser.add_argument('-i','--input', help='File or directory containing spectra in mgf/mzxml format.',required=True)
	parser.add_argument('-m','--mode', help="GNPS mode - 'library_search', 'networking', or 'both'.",required=True)
	args = parser.parse_args()
	# check if directory
	if not os.path.isdir(args.input):
		print "Please input path to a *directory* containing mzXML and/or mgf."
		sys.exit()
	in_dir = args.input
	mode = args.mode

	un = raw_input("GNPS username: ")
	pw = getpass.getpass('Password: ')
	base_url = "gnps.ucsd.edu"
	print "Logging in as user %s." %(un)

	# open output file
	out_file = open("../data/gnps.log", "w")
	out_file.write("file_name\tworkflow\tstatus\tURL\n")

	for sample in glob.glob((in_dir+"/*[mzXML|mgf]")):
		
		# get file name
		basename = sample.split("/")[(len(sample.split("/"))-1)]
		print "Processing " + basename
		ftp("upload", sample, basename, un, pw)
		
		if mode == "library_search" or mode == "both":
			# run library search
			print "Starting GNPS library search..."
			params1 = launch_workflow(un, basename, "search")
			task_id1 = invoke_workflow(base_url, params1, un, pw)
			json_results1 = wait_for_workflow_finish("gnps.ucsd.edu", task_id1)
			print "Library search %s." %(json_results1)
			# write results
			url1 = 'https://' + base_url + '/ProteoSAFe/status.jsp?task=' + task_id1
			out_file.write(basename + "\t" + "SEARCH\t"+json_results1+'\t' + url1 + '\n')
		
		if mode == "networking" or mode == "both":
			# run metabolomic networking
			print "Starting GNPS metabolomic networking..."
			params2 = launch_workflow(un, basename, "network")
			task_id2 = invoke_workflow(base_url, params2, un, pw)
			json_results2 = wait_for_workflow_finish("gnps.ucsd.edu", task_id2)
			print "Networking %s." %(json_results2)
			# write results
			url2 = 'https://' + base_url + '/ProteoSAFe/status.jsp?task=' + task_id2
			out_file.write(basename + "\t" + "NETWORK\t"+json_results2 +'\t' + url2 + '\n')
	
		# remove input file
		ftp("delete", sample, basename, un, pw)

		# write out json results
		#out_file = open("../data/gnps.json", "w")
		#results_url = 'https://' + base_url + '/ProteoSAFe/result_json.jsp?task=' + task_id + '&view=view_all_annotations_DB'
		#out_file.write(json.dumps(json.loads(requests.get(results_url, verify=False).text)))
		#out_file.close()

	out_file.close()

if __name__ == '__main__':
	main()