import ftplib
import argparse
import sys

def ftp(mode, _input, _output, server, login, password):

	# open server and login
	ftp = ftplib.FTP(server)
	try:
		ftp.login(login, password)
	except:
		print "Login failed :/"
		sys.exit()
	
	# perform ftp operation on file
	if mode == "upload":
		myfile = open(_input, "rb")
		try:
			ftp.storbinary("STOR " + _output, myfile)
			print "Upload successful."
			myfile.close()
		except:
			print "File upload failed :/"
			sys.exit()
	elif mode == "retrieve":
		myfile = open(_output, "w")
		try:
			ftp.retrbinary('RETR %s' % _input, myfile.write)
			print "Retrieve successful."
			myfile.close()
		except:
			print "File retrieve failed :/"
			sys.exit()
	elif mode == "delete":
		try:
			ftp.delete(_input)
			print "Delete successful."
		except:
			print "File delete failed :/"
			sys.exit()
	else:
		print "Mode not recognized."
		sys.exit()
	
	ftp.quit()

def main():

	parser = argparse.ArgumentParser(description='Upload, retrieve, or delete a file on a server via ftp.')
	parser.add_argument('-i','--input', help='Input file path. If uploading, input is a local file. If retrieving or deleting, input is file on server.',required=True)
	parser.add_argument('-o','--output', help="Output file path. If uploading, output is file on server. If retrieving, output is local file. Not required for 'delete' function.",required=False)
	parser.add_argument('-m','--mode', help="Mode, either 'upload', 'retrieve', or 'delete'.",required=True)
	args = parser.parse_args()
	
	# get server info from user
	server = raw_input("Server URL: ")
	un = raw_input("Enter username: ")
	pw = raw_input("Enter password: ")

	print "Logging in to %s as user %s." %(server, un)
	print "Performing %s on %s..." %(args.mode, args.output)
	ftp(args.mode, args.input, args.output, server, un, pw)

if __name__ == '__main__':
	main()