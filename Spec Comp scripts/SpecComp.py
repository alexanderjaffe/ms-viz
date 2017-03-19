from __future__ import division
from pyteomics import mzxml
import sys
import math
from operator import itemgetter
import numpy as np
import csv
from itertools import izip



input_file=sys.argv[1]

''' fast calculation of cosine similarity '''
def fast_cosine(u,v): 
    denominator = math.sqrt(np.sum(u**2) * np.sum(v**2))
    if denominator>0:
        return np.dot(u,v)/denominator
    else:
        return float("inf")


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
	msI_list= {}
	for scan in scans:
		if scan['msLevel'] == '1':
			num= scan['num']
			mzs_MSI= scan['m/z array']
			intensities= scan['intensity array']
			msI_TIC= scan['totIonCurrent'] 
			msI_list[num]= {"num": num, "mzs": mzs_MSI, "TIC": msI_TIC}
		if scan['msLevel'] == '2':
			ms1_scan_num= int(scan['precursorMz'][0]['precursorScanNum'])
			base_mz = scan['precursorMz'][0]['precursorMz']
			precursor_intensity = scan['precursorMz'][0]['precursorIntensity']
			intensity_array=scan['intensity array'].tolist()

			msI_TIC= scans[ms1_scan_num]['totIonCurrent']
			ms2_TIC= scan['totIonCurrent']
			percentComposition= precursor_intensity/msI_TIC
			for x in range(0,len(intensity_array)): 
			 	intensity_array[x]= math.log(intensity_array[x]/ms2_TIC)	
			mzs = scan['m/z array']
			num = int(scan['num'])
			percentComposition= ms2_TIC/msI_TIC
			base_peaks[num] = {"num":num, "base_mz":base_mz, "intensities":intensity_array, "mzs":mzs, "precursor_intensity": precursor_intensity,"MSI TIC": msI_TIC, "MS2 TIC": ms2_TIC}

			all_peaks = all_peaks + mzs.tolist()
	peak_min = int(math.floor(min(all_peaks)))
	peak_max = int(math.ceil(max(all_peaks)))

	#Get rid of really big variables...
	all_peaks = None
	scans = None
	r = None
	
	#Returns a list of MS2 spectra organized by scan #, and the largest and smallest precursor peaks across all MS2 in this file.
	return peak_min, peak_max, base_peaks, msI_list

""" 
	Orders Peak data by percent composition of identified compounds that have MSII 
"""
# def orderedPercentComp(new_peak_data): 
# 	percentComp_list=[]
# 	for key in new_peak_data.keys(): 
# 		percentComp_list.append((new_peak_data[key]['num'],new_peak_data[key]['Percent of Sample'],new_peak_data[key]['MS2 TIC']))
# 	lis= sorted(percentComp_list,key=lambda x: x[1], reverse=True)
# 	return lis


''' 
	F-test clustering method 
'''

	


def vectorize_peak(peak_min, peak_max, sample_data, sample_name, msI_list):	
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

	print "Running comparison..."
	# #Remove non-unique peaks; peaks that are most identical are grouped and the most intense peak from each group is kept.
	#Only compare peaks that have masses within ~3 DA of each other?
	similarities = []
	already_calculated = []
	peak_vectors_family_unique = []

	print "Finding spectra in same families..."
	f = open('sims_new', 'a+')
	for scan in peak_vectors_list:
		found = False
		#Compare to every other scan < this scan's mz + 1.5 Da
		for i in xrange(len(peak_vectors_family_unique)-1, -1, -1):
			scan2 = peak_vectors_family_unique[i]
			mass_diff = sample_data[scan]['base_mz'] - sample_data[scan2[0]]['base_mz']
			if mass_diff <= 1.5:
				#Calculate cosine similarity of these two scans' peak vectors
				sim = fast_cosine(peak_vectors[scan], peak_vectors[scan2[0]])
				f.write(str(sim) + " ")

				if sim >= 0.90:
					peak_vectors_family_unique[i].append(scan)
					found = True
					break
			else:
				break

		#Not similar to any in our list of unique peaks; add to unique list
		if not found:
			peak_vectors_family_unique.append([scan])
	f.close()
	print "done!"

	print "Finding unique peaks in sample..."
	peak_vectors_unique=[]
	g = open('sims_new_CAMS', 'a+')
	for scan in peak_vectors_list:
		found = False
		#Compare to every other scan < this scan's mz + 1.5 Da
		""" WHAT THE CORN"""
		for i in xrange(len(peak_vectors_unique)-1, -1, -1):
			scan2 = peak_vectors_unique[i]
			mass_diff = sample_data[scan]['base_mz'] - sample_data[scan2[0]]['base_mz']
			if mass_diff <= 1.5:
				#Calculate cosine similarity of these two scans' peak vectors
				weight= doCAMS(sample_data[scan], sample_data[scan2[0]])
				print weight
				g.write(str(weight) + " ")

				if weight >= 30:
					print "binned!"
					peak_vectors_unique[i].append(scan)
					found = True
					break
			else:
				break

		#Not similar to any in our list of unique peaks; add to unique list
		if not found:
			peak_vectors_unique.append([scan])
	g.close()
	print "done!"
	#Create final data for this sample, return
	#Create consensus peaks for each "compound" (group of identical scans)
	print str(len(peak_vectors_family_unique))  + " unique clustered compounds found in this sample."
	final_peaks = {}
	for scan_group in peak_vectors_family_unique:
		if len(scan_group) > 1:
			


			scan_group[] 

			consensus_peak = peak_vectors[scan_group[0]]
			msI_TIC_total= 0
			precursor_intensity_total=0

			#Get scan in this group with biggest MS1 base peak intensity
			biggest_mz = 0
			min_mz= 9999999 
			max_mz=0
			for scan in scan_group:
				temp_num= sample_data[scan]['base_mz']

				if temp_num>max_mz: 
					max_mz= temp_num
				if temp_num<min_mz: 
					min_mz=temp_num
				#msI_TIC_total+= float(sample_data[scan]['MSI TIC'])
				precursor_intensity_total+= float(sample_data[scan][precursor_intensity])
				if sample_data[scan]['base_mz'] > biggest_mz:
					biggest_mz = sample_data[scan]['base_mz']
			
			mz_num_list=[]
			ms2_intensity_list=[]
			mzs=[]
			for msI in msI_list: 
				mzs= msI['mzs']
				intensities= msI['intensity']
				for mz in mzs: 
					if (mz<=max_mz and mz=>min_mz): 
						mz_num_list.append(msI["num"])
						ms2_intensity_list.append(intensities[mzs.index(mz)])
			msI_TIC_sum=0
			ms2_intensity_sum=0 				
			for num in mz_num_list: 
				msI_TIC_sum+=msI_list[num]['TIC']
			for intensity in ms2_intensity_list: 
				ms2_intensity_sum+= intensity

			percentComp= ms2_intensity_sum/msI_TIC_sum

			




			#Create consensus spectrum
			scan1 = scan_group[0]
			scan_group.pop(0)
			for scan in scan_group:
				consensus_peak = consensus_peak + peak_vectors[scan]			

			#Re-normalize the resulting consensus spectrum
			consensus_peak = normalize(consensus_peak.reshape(1, -1), norm='l1')[0]

			peak_data = sample_data[scan1]
			peak_data['origin'] = str(peak_data['num'])
			for scan in scan_group:
				peak_data['origin'] = peak_data['origin'] + "," + str(sample_data[scan]['num'])
			peak_data['vector'] = csr_matrix(consensus_peak) #Store only as a COO matrix
			peak_data['base_mz'] = biggest_mz
			peak_data['percent composition']=percentComp
			final_peaks[scan1] = peak_data
		else:
			scan1 = scan_group[0]
			peak_data = sample_data[scan1]
			peak_data['vector'] = csr_matrix(peak_vectors[scan1]) #Store only as a COO matrix
			peak_data['origin'] = str(peak_data['num'])
			peak_data['percent composition']= peak_data['precursor_intensity']/peak_data['MSI TIC'] 


	return final_peaks




'''
Compares all ms2 between samples and creates a compound abundance table for all compounds across all samples.
'''
def compare_samples(samples_data, output_file):

	#Quit if only 1 sample...
	if len(samples_data) <= 1:
		print "There was only one sample processed! Please add more samples to your sample sheet..."
		sys.exit(1)
	#Cluster and count compounds across all samples
	compounds = []
	compounds_list = []

	#Create combined, mz sorted list of samples across all samples
	print "Creating large list..."
	for sample in samples_data:
		for scan in samples_data[sample]:
			samples_data[sample][scan]['sample'] = sample
			compounds_list.append(samples_data[sample][scan])

	compounds_list = sorted(compounds_list, key=itemgetter('base_mz'))

	for compound in compounds_list:
		found = False
		#Compare to every other scan < this scan's mz + 3 Da
		for i in xrange(len(compounds)-1, -1, -1):
			scan2 = compounds[i][0]
			mass_diff = compound['base_mz'] - scan2['base_mz']
			if mass_diff <= 1.5:
				#Calculate cosine similarity of these two scans' peak vectors
				sim = fast_cosine_csr(compound['vector'], scan2['vector'])

				if sim >= 0.90:	
					compounds[i].append(compound)
					found = True
					break
			else:
				break

		#Not similar to any in our list of unique peaks; add to unique list
		if not found:
			compounds.append([compound])

	print len(compounds)
	#Write data table to CSV

	line = "Compound,Mass"
	for sample in samples_data:
		line = line + "," + sample

	f = open(output_file, 'w+')
	f2 = open(output_file + "_spectra.txt", 'w+')

	f.write(line + "\n")
	j = 0

	for compound_group in compounds:
		masses = []
		j += 1
		for compound in compound_group:
			masses.append(compound['base_mz'])

		line = "compound_" + str(j) + "," + str(np.mean(masses)) + "+-" + str(round(np.std(masses, ddof=0),4)) + ","

		f2.write("compound_" + str(j) + "|" + str(np.mean(masses)) + "+-" + str(round(np.std(masses, ddof=0),4)))
		for compound in compound_group:
			f2.write("|" + compound['sample'] + "$" + compound['origin'])
		f2.write("\n")

		for sample in samples_data:
			found = False
			for compound in compound_group:
				if sample == compound['sample']:
					found = True
					break
			if found:
				line = line + str(compound['percent composition']) + ","
			else:
				line = line + "0,"
		line = line.rstrip(",")
		f.write(line + "\n")

	f.close()
	f2.close()



peak_min = 999999999
peak_max = 0
	#peak_data = {}
new_min, new_max, new_peak_data,msI_list = preprocess_sample(input_file)
if new_min < peak_min:
	peak_min = new_min
if new_max > peak_max:
	peak_max = new_max


peak_data = vectorize_peak(peak_min, peak_max, new_peak_data, input_file,msI_list)


#f-set

