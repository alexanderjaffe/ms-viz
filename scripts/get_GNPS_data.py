import ftplib
import argparse

def ftp(mode, _input):

	filename = "for_d3_vis.mgf"
	ftp = ftplib.FTP("ccms-ftp01.ucsd.edu")
	try:
		ftp.login("Trax_lab","Stormborn")
	except:
		print "Login failed :/"
	myfile = open(_input, "rb")
	if mode == "upload":
		try:
			ftp.storbinary("STOR " + filename, myfile)
		except:
			print "File upload failed :/"
	elif mode == "delete":
		ftp.delete(filename)
	ftp.quit()

def main():

	__author__ = "Alexander L. Jaffe"
	#parser = argparse.ArgumentParser(description='Gets chemical info from GNPS for a set of spectra.')
	#parser.add_argument('-i','--input', help='Path to directory of mzxml.',required=True)
	#args = parser.parse_args()
	#in_dir = args.input
	
	ftp("upload", "mzXML.json")
	ftp("delete", "mzXML.json")

if __name__ == '__main__':
	main()