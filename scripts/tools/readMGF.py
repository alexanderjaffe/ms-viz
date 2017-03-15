""" Reads MGF Files and returns mzs and intentsities as lists. Can definitely be modified to return more stuff"""

def readMGF(input_file): 

	mzs=[]
	intensities=[]
	with open(input_file, 'r') as f:
	        for line in f:
	        	# all m/z and intentsity values will come after the SCAN value in all mgf's 
	            if 'SCANS' in line: 
	            	for line in f:
	            		#get rid of \n variables 
	            		if line != "\n":
	            			line = line.replace("\n", "")
	            			mzs.append(float(line.split('\t')[0]))
	            			intensities.append(float(line.split('\t')[1]))

	return mzs, intensities