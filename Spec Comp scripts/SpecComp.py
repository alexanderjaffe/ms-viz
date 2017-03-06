from __future__ import division
from pyteomics import mzxml
import sys
import math
from operator import itemgetter
import numpy as np
import csv
from itertools import izip



input_file=sys.argv[1]

'''Retrieves all MS2 spectra from a sample file. log-normalize all MS2 base peak intensities.'''
def preprocess_sample(sample):
	print "Reading " + sample
	scans = []
	r = mzxml.read(sample)
	while True:
		try:
			scans.append(r.next())
		except:
			break

	print str(len(scans)) + " scans found in " + sample
	base_peaks = {}
	all_peaks = []
	for scan in scans:
		if scan['msLevel'] == '2':
			ms1_scan_num= int(scan['precursorMz'][0]['precursorScanNum'])
			base_mz = scan['precursorMz'][0]['precursorMz']
			precursor_intensity = scan['precursorMz'][0]['precursorIntensity']
			intensity_array=scan['intensity array'].tolist()

			msI_TIC= scans[ms1_scan_num]['totIonCurrent']
			ms2_TIC= scan['totIonCurrent']
			percentComposition= ms2_TIC/msI_TIC
			#Normalize and log transform peak intensities in each scan
			#intensities = normalize(np.log(1+np.asarray(scan['intensity array']).reshape(1,-1)), norm='l1')[0]
			for x in range(0,len(intensity_array)): 
			 	intensity_array[x]= intensity_array[x]/ms2_TIC	
			mzs = scan['m/z array']
			num = int(scan['num'])
			percentComposition= ms2_TIC/msI_TIC
			base_peaks[num] = {"num":num, "base_mz":base_mz, "intensities":intensity_array, "mzs":mzs, "precursor_intensity": precursor_intensity,"Percent of Sample": percentComposition, "MS2 TIC": ms2_TIC}

			all_peaks = all_peaks + mzs.tolist()

	peak_min = int(math.floor(min(all_peaks)))
	peak_max = int(math.ceil(max(all_peaks)))

	#Get rid of really big variables...
	all_peaks = None
	scans = None
	r = None
	
	#Returns a list of MS2 spectra organized by scan #, and the largest and smallest precursor peaks across all MS2 in this file.
	return peak_min, peak_max, base_peaks

""" 
	Orders Peak data by percent composition of identified compounds that have MSII 
"""
def orderedPercentComp(new_peak_data): 
	percentComp_list=[]
	for key in new_peak_data.keys(): 
		percentComp_list.append((new_peak_data[key]['num'],new_peak_data[key]['Percent of Sample'],new_peak_data[key]['MS2 TIC']))
	lis= sorted(percentComp_list,key=lambda x: x[1], reverse=True)
	return lis




''' 
	F-test clustering method 
'''

def doCAMS(peak_vector1,peak_vector2):
	fset=7
	weight=0
	for g in range(0, len(peak_vector1)-fset): 
		peak_subset1= peak_vector1[g:g+fset]
		peak_subset2= peak_vector2[g:g+fset]
		if (all(values1 == 0 for values1 in peak_subset1) and all(values2 == 0 for values2 in peak_subset2)):
			weight+=0 
		else: 	
			#to whomever has to decipher this code in the future, I apologize in advance for how awful this is
			for intensity1, intensity2 in zip(peak_subset1, peak_subset2): 
				if (intensity1==0 and intensity2==0):
					weight+=1 
				elif(intensity1==0 and intensity2!=0): 
					weight+=0 
				elif(intensity1!=0 and intensity2==0): 
					weight+=0 
				elif(abs(1-(intensity1/intensity2))<= 0.05):
					weight+=1 
				else: 
					weight+=0 
	return weight

def vectorize_peak(peak_min, peak_max, sample_data, sample_name):	
	#initializes peak vector (intervals of 0.1 so multiple all values by 10)
	vector_length = ( peak_max - peak_min ) * 10
	peak_vectors_list = []
	print "Creating peak vectors for " + sample_name

	peak_vectors = {}
	new_peaks_list = sorted(sample_data.values(), key=itemgetter('base_mz'))
	#Creates peak vectors for each scan in this sample
	for scan in new_peaks_list:
		peak_vector = np.zeros(vector_length, dtype=np.float32)
	 	i = 0
		for p in scan['mzs']:
			pos = int((math.floor(p*10)/10 - peak_min) * 10)
			peak_vector[pos] = scan['intensities'][i]
			i += 1
		peak_vectors[scan['num']] = peak_vector

	  	peak_vectors_list.append(scan['num'])
		new_peaks_list = None

	print "Finding unique peaks in sample..."
	# #Remove non-unique peaks; peaks that are most identical are grouped and the most intense peak from each group is kept.
	#Only compare peaks that have masses within ~3 DA of each other?
	similarities = []
	already_calculated = []
	peak_vectors_unique = []

	f = open('sims_new', 'a+')
	for scan in peak_vectors_list:
		found = False
		#Compare to every other scan < this scan's mz + 1.5 Da
		""" WHAT THE CORN"""
		for i in xrange(len(peak_vectors_unique)-1, -1, -1):
			scan2 = peak_vectors_unique[i]
			mass_diff = sample_data[scan]['base_mz'] - sample_data[scan2[0]]['base_mz']
			if mass_diff <= 1.5:
				#Calculate cosine similarity of these two scans' peak vectors
				weight = doCAMS(peak_vectors[scan], peak_vectors[scan2[0]])
				f.write(str(weight) + " ")

				if weight>= 30:
					peak_vectors_unique[i].append(scan)
					found = True
					break
			else:
				break

		#Not similar to any in our list of unique peaks; add to unique list
		if not found:
			peak_vectors_unique.append([scan])
	f.close()
	# #Create final data for this sample, return
	# #Create consensus peaks for each "compound" (group of identical scans)
	# print str(len(peak_vectors_unique))  + " unique clustered compounds found in this sample."
	# final_peaks = {}
	# for scan_group in peak_vectors_unique:
	# 	if len(scan_group) > 1:
	# 		consensus_peak = peak_vectors[scan_group[0]]

	# 		#Get scan in this group with biggest MS1 base peak intensity
	# 		biggest_mz = 0
	# 		for scan in scan_group:
	# 			if sample_data[scan]['base_mz'] > biggest_mz:
	# 				biggest_mz = sample_data[scan]['base_mz']

	# 		#Create consensus spectrum
	# 		scan1 = scan_group[0]
	# 		scan_group.pop(0)
	# 		for scan in scan_group:
	# 			consensus_peak = consensus_peak + peak_vectors[scan]			

	# 		#Re-normalize the resulting consensus spectrum
	# 		consensus_peak = normalize(consensus_peak.reshape(1, -1), norm='l1')[0]

	# 		peak_data = sample_data[scan1]
	# 		peak_data['origin'] = str(peak_data['num'])
	# 		for scan in scan_group:
	# 			peak_data['origin'] = peak_data['origin'] + "," + str(sample_data[scan]['num'])
	# 		peak_data['vector'] = csr_matrix(consensus_peak) #Store only as a COO matrix
	# 		peak_data['base_mz'] = biggest_mz

	# 		final_peaks[scan1] = peak_data
	# 	else:
	# 		scan1 = scan_group[0]
	# 		peak_data = sample_data[scan1]
	# 		peak_data['vector'] = csr_matrix(peak_vectors[scan1]) #Store only as a COO matrix
	# 		peak_data['origin'] = str(peak_data['num'])
	# 		final_peaks[scan1] = peak_data

	# return final_peaks


peak_min = 999999999
peak_max = 0
	#peak_data = {}
new_min, new_max, new_peak_data = preprocess_sample(input_file)
if new_min < peak_min:
	peak_min = new_min
if new_max > peak_max:
	peak_max = new_max


peak_data = vectorize_peak(peak_min, peak_max, new_peak_data, input_file)


#f-set

