import ftplib
import argparse
import sys
import requests
import json
import time
from requests.packages.urllib3.exceptions import InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

def ftp(mode, _input, login, password):

	filename = "for_d3_vis.mgf"
	ftp = ftplib.FTP("ccms-ftp01.ucsd.edu")
	try:
		ftp.login(login, password)
	except:
		print "Login failed :/"
		sys.exit()
	if mode == "upload":
		myfile = open(_input, "rb")
		try:
			ftp.storbinary("STOR " + filename, myfile)
		except:
			print "File upload failed :/"
			sys.exit()
		myfile.close()
	elif mode == "delete":
		ftp.delete(filename)
	ftp.quit()

# invokes GNPS workflow
# BY M. WANG - miw023@cs.ucsd.edu
def invoke_workflow(base_url, parameters, login, password):
	username = login
	password = password

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
		print("Launched Task: : " + r.text)
		return task_id
	else:
		print(task_id)
		return None

# Waits for not running, then returns status
# BY M. WANG - miw023@cs.ucsd.edu
def wait_for_workflow_finish(base_url, task_id):
	url = 'https://' + base_url + '/ProteoSAFe/status_json.jsp?task=' + task_id
	json_obj = json.loads(requests.get(url, verify=False).text)
	while (json_obj["status"] != "FAILED" and json_obj["status"] != "DONE"):
		#print("Waiting for task: " + task_id)
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

# creates parameter mapping for GNPS workflow
# BY M. WANG - miw023@cs.ucsd.edu
def launch_workflow():
	
	parameters_map = {}
	parameters_map["workflow"] = "MOLECULAR-LIBRARYSEARCH"
	parameters_map["email"] = "alexander_jaffe@berkeley.edu"
	parameters_map["uuid"] = "1DAB2BC6-6827-0001-9AB5-390F1E781419"
	#parameters_map["protocol"] = "none"
	parameters_map["library_on_server"] = "d.speclibs"
	#parameters_map["desc"] = "Gold"
	parameters_map["SEARCH_LIBQUALITY"] = 3
	parameters_map["spec_on_server"] = "f.ajaffe/for_d3_vis.mgf"
	parameters_map["tolerance.PM_tolerance"] = 1.5
	parameters_map["tolerance.Ion_tolerance"] = 0.5
	parameters_map["MIN_MATCHED_PEAKS"] = 6
	parameters_map["SCORE_THRESHOLD"] = 0.7
	parameters_map["ANALOG_SEARCH"] = 0
	parameters_map["MAX_SHIFT_MASS"] = 100.0
	parameters_map["TOP_K_RESULTS"] = 1
	parameters_map["FILTER_STDDEV_PEAK_INT"] = 0
	parameters_map["MIN_PEAK_INT"] = 50.0
	parameters_map["FILTER_PRECURSOR_WINDOW"] = 1
	parameters_map["FILTER_LIBRARY"] = 0
	parameters_map["WINDOW_FILTER"] = 1
	#parameters_map["TASK_ID_COMPARISON"] = task_comparison
	#parameters_map["USER_ID_COMPARISON"] = "continuous"

	return(parameters_map)

def main():

	__author__ = "Alexander L. Jaffe"
	parser = argparse.ArgumentParser(description='Gets chemical info from GNPS for a set of spectra.')
	#parser.add_argument('-i','--input', help='Path to directory of mzxml.',required=True)
	#args = parser.parse_args()
	#in_dir = args.input
	
	#un = raw_input("GNPS username: ")
	#pw = raw_input("GNPS password: ")
	un = "ajaffe"
	pw = "RivierA8"
	base_url = "gnps.ucsd.edu"

	print "Logging in as user %s." %(un)
	print "Uploading spectra file to FTP..."
	ftp("upload", "../data/test.mgf", un, pw)
	print "Invoking GNPS workflow..."
	task_id = invoke_workflow(base_url, launch_workflow(), un, pw)
	print "Submitted to GNPS with task ID %s." %(task_id)
	json_results = wait_for_workflow_finish("gnps.ucsd.edu", task_id)
	print "Task %s" %(json_results)
	ftp("delete", "for_d3_vis.mgf", un, pw)
	# write out json results
	out_file = open("../data/gnps.json", "w")
	results_url = 'https://' + base_url + '/ProteoSAFe/result_json.jsp?task=' + task_id + '&view=view_all_annotations_DB'
	out_file.write(json.dumps(json.loads(requests.get(results_url, verify=False).text)))
	#out_file.write(json.dumps(json_results))
	out_file.close()

if __name__ == '__main__':
	main()