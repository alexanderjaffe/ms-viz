from __future__ import division
from pyteomics import mgf , auxiliary
import sys
import numpy as np
import os
from operator import truediv, sub

input_file1= sys.argv[1]
input_file2= sys.argv[2]

def readMGF(input_file): 

	mzs=[]
	intensities=[]
	with open(input_file, 'r') as f:
	        for line in f:
	            if 'SCANS' in line: 
	            	for line in f:
	            		if line != "\n":
	            			line = line.replace("\n", "")
	            			mzs.append(float(line.split('\t')[0]))
	            			intensities.append(float(line.split('\t')[1]))

	return mzs, intensities




def doCAMS(input_file1,input_file2):
	fset=7
	
	weight=0 
	mzs1, intensities1=readMGF(input_file1)
	mzs2, intensities2=readMGF(input_file2)
	count=0 

	 
	for g in range(0, len(mzs1)-fset):  
		count +=1 

	 	mz_list1= mzs1[g:g+fset]
	 	intensity_list1= intensities1[g:g+fset]
	 	for h in range(0,len(mzs2)-fset): 
	 		mz_list2= mzs2[h:h+fset]
	 		intensity_list2= intensities1[h:h+fset]

	 		mass_diff= map(sub,mz_list1,mz_list2)
	 		print mass_diff
	 		raw_input("Press Enter to continue...")
	 		if (all(abs(i) <=0.3 for i in mass_diff)):
	 			intenstiy_diff= map(truediv,intensity_list1,intensity_list2)
	 			if(all(abs(1-j)<=0.05 for j in intenstiy_diff)):
	 				#print intenstiy_diff
	 				weight +=1 
	 	  
	return weight,count

print doCAMS(input_file1,input_file2)