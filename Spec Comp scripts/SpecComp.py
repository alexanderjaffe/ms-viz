from __future__ import division
from multiprocessing import Pool  
from pyteomics import mzxml
import math
import sys
import numpy as np
from scipy.sparse import csr_matrix, find
from sklearn.preprocessing import normalize
from sys import getsizeof
from operator import itemgetter
import argparse
import time 


''' fast calculation of cosine similarity (csr)'''
def fast_cosine_csr(u,v): 
    uData, vData = u.data, v.data
    denominator = math.sqrt(np.sum(uData**2) * np.sum(vData**2))
    if denominator>0:
        uCol, vCol = u.indices, v.indices 
        uI = uData[np.in1d(uCol, vCol)]
        vI = vData[np.in1d(vCol, uCol)]
        return np.dot(uI,vI)/denominator
    else:
        return float("inf")

''' fast calculation of cosine similarity '''
def fast_cosine(u,v): 
    denominator = math.sqrt(np.sum(u**2) * np.sum(v**2))
    if denominator>0:
        return np.dot(u,v)/denominator
    else:
        return float("inf")

'''Retrieves all MS2 spectra from a sample file. log-normalize all MS2 base peak intensities.'''
def preprocess_sample(sample):
	print( "Reading " + sample)
	scans = []
	r = mzxml.read(sample)
	while True:
		try:
			scans.append(r.next())
		except:
			break
	
	print( str(len(scans)) + " scans found in " + sample)
	
	if len(scans) < 200:
		print( "Taking %s sample out" %(sample))
		raise
	else: 
		print("Screening done")
	
	all_peaks = []
	
	# get frequency of ion ms2 by mz
	mz_list = []
	for scan in scans:
		if scan['msLevel'] == 2:
			base_mz = scan['precursorMz'][0]['precursorMz']
			mz_list.append(int(base_mz)*10000)

	mz_dict = {x:mz_list.count(x) for x in mz_list}

	#print( "Scanning scans in scan... "
	init_scans = 0
	end_scans = 0
	msI_list={}
	base_peaks = {}

	for scan in scans:
		if scan['msLevel'] == 1:
			num= int(scan['num'])
			intensity_array_MSI= scan['intensity array']
			mzs_MSI= scan['m/z array']
			TIC_MSI= scan['totIonCurrent']
			msI_list[num]= {"num": num,"TIC": TIC_MSI, "mzs": mzs_MSI, "intensity": intensity_array_MSI}
		if scan['msLevel'] == 2:
			
			base_mz = scan['precursorMz'][0]['precursorMz']
			init_scans +=1
			#print mz_dict[(int(base_mz)*10000)] in mz_dict
			#filter mzs to keep ions w/ > 10 occurrences
			if mz_dict[(int(base_mz)*10000)] > 10:				
				
				end_scans +=1
				num2 = int(scan['num'])
				ms1_scan_num= int(scan['precursorMz'][0]['precursorScanNum'])
				msI_TIC= scans[ms1_scan_num]['totIonCurrent']	
				precursor_intensity= scan['precursorMz'][0]['precursorIntensity']
				mzs_MS2 = scan['m/z array']
				intensity_array_ms2_temp= scan['intensity array']
				intensity_array_ms2= intensity_array_ms2_temp/float(msI_TIC)

				base_peaks[num2] = {"num":num2, "base_mz":base_mz, "intensities":intensity_array_ms2, "mzs":mzs_MS2, "MSI TIC": msI_TIC,"precursor_intensity": precursor_intensity}
				all_peaks = all_peaks + mzs_MS2.tolist()

	peak_min = int(math.floor(min(all_peaks)))
	peak_max = int(math.ceil(max(all_peaks)))

	#"Get rid of really big variables and hoping to the MS Gods that this doesn't shut down..." 
	all_peaks = None
	scans = None
	r = None
	mz_list=None

	print( "%d scans filtered out by ms2 freq" %(init_scans - end_scans))
	#Returns a list of MS2 spectra organized by scan #, and the largest and smallest precursor peaks across all MS2 in this file.
	return peak_min, peak_max, base_peaks,msI_list



'''
(1) Converts all of the MS2 peak lists in a sample to a vectorized format.
(2) Removes non-unique MS2 spectra (precursor mass delta of 3 mz; cosine similarity > 0.97)
(3) Creates a weighted consensus spectra for near-identical MS2 spectra
'''
def vectorize_peak(peak_min, peak_max, sample_data, sample_name, msI_list):	

	#initializes peak vector (intervals of 0.1 so multiple all values by 10)
	vector_length = ( peak_max - peak_min ) * 10
	peak_vectors_list = []
	print( "Creating peak vectors for " + sample_name)

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

	print( "Finding unique peaks in sample...")
	# #Remove non-unique peaks; peaks that are most identical are grouped and the most intense peak from each group is kept.
	#Only compare peaks that have masses within ~3 DA of each other?
	similarities = []
	already_calculated = []
	peak_vectors_unique = []


	f = open('sims_new', 'a+')
	for scan in peak_vectors_list:
		found = False
		#Compare to every other scan < this scan's mz + 1.5 Da
		for i in xrange(len(peak_vectors_unique)-1, -1, -1):
			scan2 = peak_vectors_unique[i]
			
			mass_diff = sample_data[scan]['base_mz'] - sample_data[scan2[0]]['base_mz']
			
			if mass_diff <= 1.5:
				#Calculate cosine similarity of these two scans' peak vectors
				sim = fast_cosine(peak_vectors[scan], peak_vectors[scan2[0]])
				f.write(str(sim) + " ")

				if sim >= 0.90:
					peak_vectors_unique[i].append(scan)
					found = True
					break
			else:
				break

		#Not similar to any in our list of unique peaks; add to unique list
		if not found:
			peak_vectors_unique.append([scan])
	f.close()

	#Create final data for this sample, return
	#Create consensus peaks for each "compound" (group of identical scans)
	print( str(len(peak_vectors_unique))  + " unique clustered compounds found in this sample.")
	final_peaks = {}
	print("Creating consensus peaks...")
	count=0 
	for scan_group in peak_vectors_unique:
		if len(scan_group) > 1:
			consensus_peak = peak_vectors[scan_group[0]]

			#Get scan in this group with biggest MS1 base peak intensity
			biggest_mz = 0
			min_mz= 9999999 
			max_mz=0

			#creates range of mz values within group to compare MSI mzs to 
			for scan in scan_group:
				temp_num= sample_data[scan]['base_mz']
				if temp_num>max_mz: 
					max_mz= temp_num
				if temp_num<min_mz: 
					min_mz=temp_num
				if sample_data[scan]['base_mz'] > biggest_mz:
					biggest_mz = sample_data[scan]['base_mz']
			
			max_mz+= 0.00005
			min_mz-= 0.00005 		

			mz_num_list=[]
			ms2_intensity_list=[]





			for msI, msI_items in msI_list.iteritems(): 	
				for mz in msI_items["mzs"]: 
					if (mz<=max_mz and mz>=min_mz): 
						mz_num_list.append(msI_items["TIC"])
						index=  np.extract(msI_items["mzs"]==mz, msI_items["mzs"])
						ms2_intensity_list.append(msI_items["intensity"][index[0]])
					elif(mz> max_mz):
						break
			if len(mz_num_list)==0: 
				print(max_mz, min_mz)
				for scan in scan_group: 
					print (sample_data[scan]["base_mz"])	
			msI_TIC_sum=0
			ms2_intensity_sum=0 			
			#print mz_num_list,scan_group

			msI_TIC_sum=sum(mz_num_list) 
			for intensity in ms2_intensity_list: 
				ms2_intensity_sum+= intensity
			
			percentComp= ms2_intensity_sum/float(msI_TIC_sum)
			

			#Create consensus spectrum
			scan1 = scan_group[0]
			scan_group.pop(0)
			for scan in scan_group:
				consensus_peak = consensus_peak + peak_vectors[scan]			


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
			final_peaks[scan1] = peak_data
		count +=1
		print "Scan group clustered: " + str(count)
	return final_peaks

'''
Compares all ms2 between samples and creates a compound abundance table for all compounds across all samples.
'''
def compare_samples(samples_data, output_file):

	#Quit if only 1 sample...
	if len(samples_data) <= 1:
		print( "There was only one sample processed! Please add more samples to your sample sheet...")
		sys.exit()
	#Cluster and count compounds across all samples
	compounds = []
	compounds_list = []

	#Create combined, mz sorted list of samples across all samples
	print( "Creating large list...")
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

	print( len(compounds))
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

def main():


	#Read in sample data from mapping file
	__author__ = "Alex Crits-Christoph, Alexander Jaffe AND RYAN NGUYEN SUCK IT,MA HUNPS"
	parser = argparse.ArgumentParser(description='Processes a list of mzXML files as described in a mapping file to cluster and compare spectra across samples.')
	parser.add_argument('-i','--input', help='Path to input mapping file (tab separated, "file	sample	grouping" 3+ column header',required=True)
	parser.add_argument('-f','--filter_singletons', help="Do not include compounds that aren't replicated by at least 2 scans for a given sample. (similar to GNPS clusters)",required=False)

	parser.add_argument('-o','--output', help='Name of output file (default: compound_table.txt)',required=False)

	args = parser.parse_args()

	if not args.output:
		output_file = "compound_table.txt"
	else:
		output_file = args.output

	mapping_file = args.input

	#Get the list of sample files from the mapping file.
	samples = []
	f = open(mapping_file)
	i = 0 
	for line in f.readlines():
		if i != 0:
			samples.append(line.split("\t")[0])
		else:
			i += 1
	f.close()

	#Preprocess all samples

	#Get min and max peaks across all MS2 in all samples in our dataset.
	peak_min = 999999999
	peak_max = 0
	peak_data = {}
	msI_list={}

	#Get MS2 peak data for all samples in this dataset.
	print("Now Getting MS2...")
	for sample in samples:
		try:
			new_min, new_max, new_peak_data,msI_data = preprocess_sample(sample)
		except:
			continue

		if new_min < peak_min:
			peak_min = new_min
		if new_max > peak_max:
			peak_max = new_max
		peak_data[sample] = new_peak_data
		msI_list[sample]= msI_data  

		new_peak_data=None 
		msI_data=None  

	for sample in peak_data:
		peak_data[sample] = vectorize_peak(peak_min, peak_max, peak_data[sample], sample,msI_list[sample])

	#Compare samples
	print( "Comparing samples..."	)
	compare_samples(peak_data, output_file)

if __name__ == '__main__':
	main()
