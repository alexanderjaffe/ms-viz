import ftplib
import argparse
import sys

def ftp(mode, _input, login, password):

	filename = "for_d3_vis.mgf"
	ftp = ftplib.FTP("ccms-ftp01.ucsd.edu")
	try:
		ftp.login(login, password)
	except:
		print "Login failed :/"
		sys.exit()
	myfile = open(_input, "rb")
	if mode == "upload":
		try:
			ftp.storbinary("STOR " + filename, myfile)
		except:
			print "File upload failed :/"
			sys.exit()
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
        'user' : username,
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
        print("Waiting for task: " + task_id)
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
    parameters_map["protocol"] = "none"
    parameters_map["library_on_server"] = "d.speclibs"
    parameters_map["desc"] = desc_string
    parameters_map["spec_on_server"] = spec_string
    parameters_map["tolerance.PM_tolerance"] = pm_tolerance
    parameters_map["tolerance.Ion_tolerance"] = fragment_tolerance
    parameters_map["MIN_MATCHED_PEAKS"] = min_matched_peaks
    parameters_map["SCORE_THRESHOLD"] = score_threshold
    parameters_map["ANALOG_SEARCH"] = analog_search
    parameters_map["MAX_SHIFT_MASS"] = max_mass_shift
    parameters_map["TOP_K_RESULTS"] = top_k_results
    parameters_map["FILTER_STDDEV_PEAK_INT"] = std_dev_filter
    parameters_map["MIN_PEAK_INT"] = min_peak_int
    parameters_map["FILTER_PRECURSOR_WINDOW"] = filter_precursor_window
    parameters_map["FILTER_LIBRARY"] = filter_library
    parameters_map["WINDOW_FILTER"] = window_filter
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

	print "Logging in as user %s." %(un)
	print "Uploading spectra file to FTP..."
	ftp("upload", "mzXML.json", un, pw)
	print "Invoking GNPS workflow..."
	parameter_map = launch_workflow()
	#task_id = invoke_workflow("gnps.ucsd.edu", parameters_map, username, password)
	ftp("delete", "mzXML.json", un, pw)

if __name__ == '__main__':
	main()