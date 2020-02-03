''' Function that records which binding site a lipid is interacting at, and how many contacts it makes when there

There is functionality to add s.t.  one can measure whether beads in contact are charged (positive or negative), polar, apolar or non-polar'''

import sys
import re
import ast
import numpy
import time
import MDAnalysis
from MDAnalysis.core.parallel.distances import distance_array
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.ticker as ticker

import ast
import argparse

from interaction_freq_core import protein_info, get_real_resnums_strlist, bondini_residuesymbol_size, martini_residuesymbol_size

lipid_particles = {
	'headgroup' : {'PIP2': ['PO3','PO1','PO2','RP1','RP2','RP3'],
			'PI3': ['PO3','PO0','PO1','PO2','RP1','RP2','RP3'],
			'PSPI' : ['PO3','RP1','RP2','RP3'],
			'CHOL' : ['ROH', 'R1', 'R2', 'R3', 'R4', 'R5'], #'C1', 'C2'], # this is actually the whole CHOL},
			'GM3': ['INA','B1A','B2A','B3A','B1B','B2B','B3B','INB','B1C','B2C','B3C','INC','B4C','B5C'],
			'PPCS' : ['PO4','NC3'],
			'POPS' : ['PO4','CNO'],
			'POPC' : ['PO4','NC3'],
			'POPE' : ['PO4','NH3']},
	'phosphate' : {'PIP2': ['PO3','PO1','PO2'],
			'PI3': ['PO3','PO0','PO1','PO2'],
			'PSPI' : ['PO3'],
			'CHOL' : ['ROH'],
			'GM3': ['B1A','B2A','B3A','INA'], # This is obviously not a phosphate group, but it's in the position of one in the memb
			'PPCS' : ['PO4'],
			'POPS' : ['PO4'],
			'POPC' : ['PO4'],
			'POPE' : ['PO4']},
	'phosphate_singlebead' : {'PIP2': ['PO3'],
			'PI3': ['PO3'],
			'PSPI' : ['PO3'],
			'CHOL' : ['ROH'],
			'GM3': ['INA'], # This is obviously not a phosphate group, but it's in the position of one in the memb
			'PPCS' : ['PO4'],
			'POPS' : ['PO4'],
			'POPC' : ['PO4'],
			'POPE' : ['PO4']}}
edge_colors = ['orange','blue','darkgrey']
sys_colors = {'PM':'orange', \
              'PM_whole':'orange', \
              'PEnP2':'royalblue', \
              'PC':'grey', \
              'PC_noKir':'violet', \
              'noCHOL' : 'cadetblue', \
              'PM_3x3':'orange', \
              'noGM3P2_3x3':'turquoise', \
              'PM_6x6':'y', \
              'PC_3x3':'grey', \
              'noP2_3x3' : 'royalblue', \
              'noGM3_3x3' : 'green', \
              'noCHOL_3x3' : 'mediumorchid', \
              'noP2GM3CHOL_3x3' : 'thistle', \
              'noPS_3x3':'red', \
              'PIP3_3x3' :'red', \
              'PM_noprot':'red'}
files_dict = {
	'PM':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPI_PM_3x3_ARCHER/processed_trajectories/system_3x3_fromem1x1.noW.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPI_PM_3x3_ARCHER/processed_trajectories/md_prod3x3.until20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPI_PM_3x3_ARCHER/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'bs_info_path':'/sansom/n42/bioc1280/KIR_lipids_analysis/correlation_analysis/PM/', 'nprot':9, 'nframes':50000, 'dt':1},
	'P2subPE':{
		'GRO':'/sansom/n42/bioc1280/Kir2_2_3SPI_PIP2subPE_3x3/processed_trajectories/md_prod3x3.noW.firstframe.gro',
		'XTC':'/sansom/n42/bioc1280/Kir2_2_3SPI_PIP2subPE_3x3/processed_trajectories/md_prod3x3.until10us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir2_2_3SPI_PIP2subPE_3x3/processed_trajectories/md_prod3x3.from10us.to20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir2_2_3SPI_PIP2subPE_3x3/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'bs_info_path':'/sansom/n42/bioc1280/KIR_lipids_analysis/correlation_analysis/P2subPE/', 'nprot':9, 'nframes':50000, 'dt':1},
	'PMnoGM3':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noGm3_3x3/processed_trajectories/md_prod3x3.noW.firstframe.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noGm3_3x3/processed_trajectories/md_prod3x3.until10us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPIsys3_noGm3_3x3/processed_trajectories//md_prod3x3.from10us.to20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPIsys3_noGm3_3x3/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'bs_info_path':'/sansom/n42/bioc1280/KIR_lipids_analysis/correlation_analysis/PMnoGM3/', 'nprot':9, 'nframes':50000, 'dt':1},
	'PMnoCHOL':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noCHOL3_3x3/processed_trajectories/md_prod3x3.until20us.noW.firstframe.pbcmol.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noCHOL3_3x3/processed_trajectories/md_prod3x3.until20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPIsys3_noCHOL3_3x3/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'nprot':9, 'nframes':50000, 'dt':1},
	'PMnoP2GM':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPI_noP2noGM3_3x3_ARCHER/processed_trajectories/md_prod3x3.until20us.noW.firstframe.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPI_noP2noGM3_3x3_ARCHER/processed_trajectories/md_prod3x3.until20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPI_noP2noGM3_3x3_ARCHER/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'bs_info_path':'/sansom/n42/bioc1280/KIR_lipids_analysis/correlation_analysis/PMnoP2GM/', 'nprot':9, 'nframes':50000, 'dt':1},
	'PMnoP2GMCh':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noP2noGm3noCHOL2_3x3/processed_trajectories/md_prod3x3.until20us.noW.firstframe.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noP2noGm3noCHOL2_3x3/processed_trajectories/md_prod3x3.until20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPIsys3_noP2noGm3noCHOL2_3x3/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'nprot':9, 'nframes':50000, 'dt':1},
	'P2noPS':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPIsys3_P2noPS_3x3/processed_trajectories/md_prod3x3_P2noPS.from0us.to20us.noW.firstframe.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPIsys3_P2noPS_3x3/processed_trajectories/md_prod3x3_P2noPS.from0us.to20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPIsys3_P2noPS_3x3/processed_trajectories/md_prod3x3_P2noPS.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'bs_info_path':'/sansom/n42/bioc1280/KIR_lipids_analysis/correlation_analysis/P2noPS/', 'nprot':9, 'nframes':50000, 'dt':1},
	'PM_PIP3':{
		'GRO':'/sansom/n16/bioc1280/curie_prodruns_3SPI_3x3_PIP3/processed_trajectories/md_prod3x3.noW.firstframe.gro',
		'XTC':'/sansom/n16/bioc1280/curie_prodruns_3SPI_3x3_PIP3/processed_trajectories/md_prod3x3.noW.dt1ns.pbcmol.center.xtc',
		'nprot':9, 'nframes':50000, 'dt':1},
		
	'sysname':{
		'GRO':'',
		'XTC':'',
		'nprot':0, 'nframes':0, 'dt':0},
	'':{
		'GRO':'',
		'XTC':'',
		'nprot':0, 'nframes':0, 'dt':0}
	}
	

def find_prot_residues(res_total, residue_list, nrepeats):
    '''res_total is the total number of residues in one repeat, residue list should be integers
     this function will find res IDs of all the repeats in the system, 
     and write selection string for them'''
    residue_dictionary = dict( zip(residue_list, [[] for i in residue_list]) )
    for residue in residue_list:
        for repeat in range(nrepeats):
            residue_dictionary[residue].append(int(residue) + repeat*res_total - 1) # the minus 1 is because residue indexes start at 0
    return residue_dictionary

def make_split_list_single_prot(prot_seq, ff_ressize_dict=martini_residuesymbol_size):
    splitlist_single = []
    counter = 0
    for res in prot_seq:
        counter += ff_ressize_dict[res]
        splitlist_single.append(counter)
    return splitlist_single[:-1]  

def make_res_list(res_string, prot_seq):
	# residue list should be input in the format 1:4,6:9,300,310:333
	# the following decodes that format into a straight list
	if res_string == 'all':
		res_list = numpy.arange(len(list(prot_seq))) + 1
	else:
		res_list_ranges = [x.split(':') for x in res_string.split(',')]
		res_list = []
		for res_range in res_list_ranges:
			res_list += range(int(res_range[0]), int(res_range[-1])+1 ) 
	#print 'res_list: ', res_list
	return res_list

def count_contacts(gro_file, trr_file, prot_seq, protein_residue_list, nprots, lipid, bs_dict, file_label, nframes, stride=1, lipid_part='headgroup', cutoff=6.5):
	universe = MDAnalysis.Universe(gro_file, trr_file)

	protein_res_total = len(prot_seq)

	lipid_selection = 'resname {} and ('.format(lipid)
	for bead in lipid_particles[lipid_part][lipid]:
		lipid_selection += 'name {} or '.format(bead)
	lipid_selection = lipid_selection[:-4] + ')'

	lipids = universe.selectAtoms(lipid_selection)
	n_lipid_beads = len(lipid_particles[lipid_part][lipid])
	n_lipids = lipids.numberOfAtoms() / n_lipid_beads

    #initialise protein-lipid interactions frequency list
	# initialise data storage
	n_contacts_dict = {}
	n_contacts_dict_pertime_avg = {}
	n_contacts_dict_pertime_stdev = {}
	for bs_annotation in bs_dict.keys():
		n_contacts_dict[bs_annotation] = []
		n_contacts_dict_pertime_avg[bs_annotation] = {}
		n_contacts_dict_pertime_stdev[bs_annotation] = {}
		for prot in range(nprots):
			n_contacts_dict_pertime_avg[bs_annotation][prot] = []
			n_contacts_dict_pertime_stdev[bs_annotation][prot] = []
	#print bs_dict, n_contacts_dict
	
	startTime = time.time()    
	print 'Here we go...'
	frame = 0
	for ts in universe.trajectory[:nframes+1:stride]:
		if frame >= nframes:
			print 'Have reached maximum number of frames specified.  Stopping...'
			continue
		for i in range(nprots):
			single_prot = universe.segments[0][range(i*protein_res_total,(i+1)*protein_res_total)]
			for bs_annotation in bs_dict.keys():
				#print 'i, bs_annotation', i, bs_annotation
				bs_residues = bs_dict[bs_annotation]
				repeat_res_list = i*protein_res_total + numpy.array(bs_residues)
				single_prot_bs_coords = universe.segments[0][repeat_res_list].coordinates()
				all_dists = distance_array(single_prot_bs_coords, lipids.coordinates(), ts.dimensions)
				#print 'all_dists.shape', all_dists.shape
				#print 'bs_residues, prot_seq[bs_residues]', bs_residues, numpy.array(list(prot_seq))[bs_residues]
				prot_split_bs = make_split_list_single_prot(numpy.array(list(prot_seq))[bs_residues])
				protein_lipid_dist_perresidue_all = numpy.array([[x.min() for x in numpy.split(lip, prot_split_bs, axis=0)] for lip in numpy.split(all_dists, n_lipids, axis=1)])
				#print 'len(bs_residues), protein_lipid_dist_perresidue_all.shape', len(bs_residues), protein_lipid_dist_perresidue_all.shape
				bs_interactions = protein_lipid_dist_perresidue_all <= cutoff
				bs_interactions_ncontacts_perlipid = numpy.sum(bs_interactions, axis = 1)
				#print 'bs_interactions_ncontacts_perlipid.shape', bs_interactions_ncontacts_perlipid.shape
				if numpy.sum(bs_interactions_ncontacts_perlipid) > 0:
					in_contact_bs_interactions_ncontacts_perlipid = bs_interactions_ncontacts_perlipid[bs_interactions_ncontacts_perlipid>0]
					#print 't, pr, bs, all_dists.sh, in_ctct.sh', frame, i, bs_annotation, all_dists.shape, in_contact_bs_interactions_ncontacts_perlipid.shape
					n_contacts_dict[bs_annotation].append(list(in_contact_bs_interactions_ncontacts_perlipid))
					n_contacts_dict_pertime_avg[bs_annotation][i].append(numpy.mean(in_contact_bs_interactions_ncontacts_perlipid))
					n_contacts_dict_pertime_stdev[bs_annotation][i].append(numpy.std(in_contact_bs_interactions_ncontacts_perlipid))
				else:
					n_contacts_dict_pertime_avg[bs_annotation][i].append(0)
					n_contacts_dict_pertime_stdev[bs_annotation][i].append(0)
				#print n_contacts_dict
		frame += 1
		if frame % 100 == 0:
			print 'Frame {} took {:3f} s\r'.format(frame, time.time()-startTime)
		startTime = time.time()
	f = open('n_contacts_{}.txt'.format(file_label), 'w')
	f.write(str(n_contacts_dict)+'\n\n')
	f.write(str(n_contacts_dict_pertime_avg)+'\n\n')
	f.write(str(n_contacts_dict_pertime_stdev)+'\n\n')
	f.close()
	return n_contacts_dict, n_contacts_dict_pertime_avg, n_contacts_dict_pertime_stdev

#function to write data from dict, avg and std

# function to read data from simulations of different systems
def read_ncontacts_data(filename):
	f = open(filename)
	n_contacts_dict = ast.literal_eval(f.readlines()[0])
	f.close()
	return n_contacts_dict

# make lists of binding sites, by index, with n residues per bs and bs colour scheme:
def parse_binding_sites(sys_path,system,lipid,fitting, annotations_list, n_monomers, nres_monomer):
	path  = sys_path + '/' + lipid + '/'
	ordered_annotations_list = []
	ordered_annotations_colours_list =[]
	f = open(path+'restimes_on_struct_restimes_and_colours_{}_{}_{}fitting.txt'.format(system, lipid, fitting))
	for line in f.readlines():
			m = re.match(r'\(\'(\d+)\', [\d\.]+\) (\w+)',line)
			if m:
				if annotations_list == []:
						ordered_annotations_list.append(int(m.group(1)))
						ordered_annotations_colours_list.append(m.group(2))
				else:
					annotation = m.group(1)
					if annotation in annotations_list:
						ordered_annotations_list.append(int(annotation))
						ordered_annotations_colours_list.append(m.group(2))
	f.close()
	print 'ordered_annotations_list, ordered_annotations_colours_list', ordered_annotations_list, ordered_annotations_colours_list
	# create dict of resIDs for each annotation in ordered_annotations_list
#	ordered_byannot_resIDs_list = []
	bs_annotation_resIDs_dict = {}
	for i,annotation in enumerate(ordered_annotations_list):
		af = open(path+'modularity_network_{}.txt'.format(annotation))
		line_of_interest = af.readlines()[1]
		reslist = numpy.array(ast.literal_eval(line_of_interest)) # + 1 # +1 is because MDanalysis reads 1st residue as residue 0, whereas gro file has 1st residue as residue 1
		whole_list = []
		for monomer in range(n_monomers):
			whole_list += list(reslist + nres_monomer*monomer)
#		ordered_byannot_resIDs_list.append(whole_list)
		bs_annotation_resIDs_dict[annotation] = whole_list
		af.close()
	#print ordered_byannot_resIDs_list
	return ordered_annotations_list, ordered_annotations_colours_list, bs_annotation_resIDs_dict

# function to do a violin plot of n contacts for each binding site
def plot_n_contacts_violin(n_contacts_dict, ordered_annotations_colours_list, ordered_annotations_list, file_label=''):
	all_n_contacts = []
	for i,bs in enumerate(ordered_annotations_list):
		n_contacts_bs = []
		for n_contacts in n_contacts_dict[bs]:
			n_contacts_bs += n_contacts 
		all_n_contacts.append(n_contacts_bs)
	fig, ax = plt.subplots(1,1, figsize=(9,4))
	#ax.boxplot(all_dists)
	#ax.boxplot(df)
	dist_violin = plt.violinplot(all_n_contacts,showmeans=True)
	#ax.violinplot(df)
	for i,pc in enumerate(dist_violin['bodies']):
		pc.set_facecolor(ordered_annotations_colours_list[i])
		pc.set_edgecolor(ordered_annotations_colours_list[i])
	plt.setp(dist_violin['cbars'], color='black')
	plt.setp(dist_violin['cmins'], color='black')
	plt.setp(dist_violin['cmaxes'], color='black')
	plt.xlabel('Binding site', fontsize=16)
	plt.ylabel('Number of residue - lipid contacts per lipid', fontsize=16)
	#plt.xticks(fontsize=16,rotation=45, ha='right')
	plt.yticks(fontsize=16)
	plt.savefig('N_contacts{}.svg'.format(file_label),bbox_inches='tight')
	plt.close
def plot_n_contacts(n_contacts_dict, ordered_annotations_colours_list, ordered_annotations_list, file_label=''):
	all_n_contacts = []
	max_n_contacts = 0
	for i,bs in enumerate(ordered_annotations_list):
		n_contacts_bs = []
		for n_contacts in n_contacts_dict[bs]:
			n_contacts_bs += n_contacts 
		all_n_contacts.append(n_contacts_bs)
		if n_contacts_bs:
			if max_n_contacts < max(n_contacts_bs):
				max_n_contacts = max(n_contacts_bs)
	fig, ax = plt.subplots(1,len(ordered_annotations_list), figsize=(2.5*len(ordered_annotations_list),2))
	for i,annotation in enumerate(ordered_annotations_list):
		if all_n_contacts[i] == []:
			continue
		else:
			print i, ordered_annotations_colours_list[i]
			barcol = ordered_annotations_colours_list[i]
			if barcol == 'green':
				barcol = 'lime'
			ax[i].bar(range(1,max_n_contacts+1),[all_n_contacts[i].count(n) for n in range(1,max_n_contacts+1)], color=barcol,align='center')
			#ax[i].set_xticks(ha='right')
			#ax[i].set_yticks(fontsize=10)
#			for tick in ax[i].get_xticklabels():
#				tick.set_fontname("Helvetica")
#			for tick in ax[i].get_yticklabels():
#				tick.set_fontname("Arial")
			ax[i].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
			ax[i].yaxis.major.formatter._useMathText = True
			ax[i].xaxis.set_major_locator(ticker.MultipleLocator(3))
			ax[i].set_xlim(0,18)
	plt.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.9, hspace=0.2, wspace=0.3)
	#ax.boxplot(all_dists)
	#ax.boxplot(df)
	#dist_violin = plt.violinplot(all_n_contacts,showmeans=True)
	#ax.violinplot(df)
#	for i,pc in enumerate(dist_violin['bodies']):
#		pc.set_facecolor(ordered_annotations_colours_list[i])
#		pc.set_edgecolor(ordered_annotations_colours_list[i])
#	plt.setp(dist_violin['cbars'], color='black')
#	plt.setp(dist_violin['cmins'], color='black')
#	plt.setp(dist_violin['cmaxes'], color='black')
#	plt.xlabel('Binding site', fontsize=16)
#	plt.ylabel('Number of residue - lipid contacts per lipid', fontsize=16)
	#def _1gaussian(x, amp1,cen1,sigma1):
	#	return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2)))
	#popt_gauss, pcov_gauss = scipy.optimize.curve_fit(_1gaussian, x_array, y_array_gauss, p0=[amp1, cen1, sigma1])
	#perr_gauss = np.sqrt(np.diag(pcov_gauss))
	#def _2gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2):
	#	return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2))) + amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen2)/sigma2)**2)))
	#popt_2gauss, pcov_2gauss = scipy.optimize.curve_fit(_2gaussian, x_array, y_array_2gauss, p0=[amp1, cen1, sigma1, amp2, cen2, sigma2])
	#perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
	#pars_1 = popt_2gauss[0:3]
	#pars_2 = popt_2gauss[3:6]
	#gauss_peak_1 = _1gaussian(x_array, *pars_1)
	#gauss_peak_2 = _1gaussian(x_array, *pars_2)
	plt.savefig('N_contacts_bar_{}.svg'.format(file_label),bbox_inches='tight')
	plt.savefig('N_contacts_bar_{}.pdf'.format(file_label),bbox_inches='tight')
	plt.close

# main
def main_get_data(gro_file, trr_file, prot_seq, protein_residue_list, nprots, lipid, file_label, sys_path, system, fitting, annotations_list, n_monomers, nres_monomer, nframes, stride=1, lipid_part='headgroup', cutoff=6.5):
	print 'Parsing binding site data...'
	ordered_annotations_list, ordered_annotations_colours_list, bs_annotation_resIDs_dict = parse_binding_sites(sys_path,system,lipid,fitting, annotations_list, n_monomers, nres_monomer)
	print 'Assessing n contacts...'
	n_contacts_dict, n_contacts_dict_pertime_avg, n_contacts_dict_pertime_stdev = count_contacts(gro_file, trr_file, prot_seq, protein_residue_list, nprots, lipid, bs_annotation_resIDs_dict, file_label, nframes, stride, lipid_part, cutoff)
	print 'Creating graph...'
	plot_n_contacts(n_contacts_dict, ordered_annotations_colours_list, ordered_annotations_list, file_label)

def main_plot_data(data_filename, sys_path,system,lipid,fitting, annotations_list, n_monomers, nres_monomer, file_label):
	print 'Parsing binding site data...'
	ordered_annotations_list, ordered_annotations_colours_list, bs_annotation_resIDs_dict = parse_binding_sites(sys_path,system,lipid,fitting, annotations_list, n_monomers, nres_monomer)
	print 'Parsing n contacts from existing data file...'
	n_contacts_dict = read_ncontacts_data(data_filename)
	print 'Creating graph...'
	plot_n_contacts(n_contacts_dict, ordered_annotations_colours_list, ordered_annotations_list, file_label)

	
# parse input
if __name__ == "__main__":
	def process_data(args):
		f = open('inputs_ncontacts_process_data_bs.log', 'w')
		f.write('Args used: {}'.format(args)+'\n'+'command used: python'+' '.join(sys.argv))
		f.close()

		# set nres_monomer and n_monomers - if numbers are given in command line this will override info in protein_info
		system = args.sys
		protname = args.p
		#if args.permonomer and protname not in monomer_format_dict:
		#	print 'Add {} to monomer_format_dict'
		#elif args.permonomer:
		#	monomer_format = monomer_format_dict[protname]
		#else: 
		#	monomer_format = None
		if protname in protein_info and args.nresm == 0:
			nres_monomer = len(protein_info[protname]['protein_seq'])
		else:
			nres_monomer = args.nresm
		if protname in protein_info and args.nresm == 0:
			n_monomers = protein_info[protname]['protein_nmonomers']
		else: 
			n_monomers = args.nmonomer
		if protname in protein_info and args.seq == '':
			prot_seq = protein_info[protname]['protein_seq']*n_monomers
		else:
			prot_seq = args.seq * n_monomers

		res_list = make_res_list(args.reslist, prot_seq)

		if args.f == 'none' and system not in files_dict:
			print 'Add gro filename to files_dict or specify using the -f flag'
		elif args.f == 'none':
			gro_file = files_dict[system]['GRO']
		else:
			gro_file = args.f
		if args.x == 'none' and system not in files_dict:
			print 'Add xtc filename(s) to files_dict or specify using the -x flag'
		elif args.x == 'none':
			xtc_file = list(files_dict[system]['XTC'].split(' '))
			nprot = files_dict[system]['nprot']
	#		dt = files_dict[system]['dt']/args.stride
			sys_path = files_dict[system]['bs_info_path']
		else:
			xtc_file = args.x
			nprot = args.np
	#		dt = args.dt
			sys_path = args.syspath
		if args.nframes == 0 and system not in files_dict:
			print 'Add nframes to files_dict or specify using the -nframes flag'
		elif args.nframes == 0:
			nframes = files_dict[system]['nframes']/args.stride
		else:
			nframes = args.nframes

		if args.fl == '':
			file_label = system + '_' + args.lipid
		else:
			file_label = args.fl

		main_get_data(gro_file, xtc_file, prot_seq, res_list, nprot, args.lipid, file_label, sys_path, system, args.fitting, args.annotations, n_monomers, nres_monomer, nframes, args.stride, args.lp, args.d)

	def plot_data(args):
		f = open('inputs_ncontacts_plot_data_bs.log', 'w')
		f.write('Args used: {}'.format(args)+'\n'+'command used: python'+' '.join(sys.argv))
		f.close()

		# set nres_monomer and n_monomers - if numbers are given in command line this will override info in protein_info
		system = args.sys
		protname = args.p
		if protname in protein_info and args.nresm == 0:
			nres_monomer = len(protein_info[protname]['protein_seq'])
		else:
			nres_monomer = args.nresm
		if protname in protein_info and args.nresm == 0:
			n_monomers = protein_info[protname]['protein_nmonomers']
		else: 
			n_monomers = args.nmonomer

		if args.syspath == 'none' and system not in files_dict:
			print 'Add syspath to files_dict or specify using the -syspath flag'
		elif args.syspath == 'none':
			sys_path = files_dict[system]['bs_info_path']
		else:
			sys_path = args.syspath

		if args.fl == '':
			file_label = system + '_' + args.lipid
		else:
			file_label = args.fl

		if args.data == '':
			data_filename = 'n_contacts_{}.txt'.format(file_label)
		else:
			data_filename = args.data

		main_plot_data(data_filename, sys_path, system, args.lipid, args.fitting, args.annotations, n_monomers, nres_monomer, file_label)

	parser = argparse.ArgumentParser(description='Run n contacts on each modularity network', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = parser.add_subparsers()

	parser_data = subparsers.add_parser('data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_data.add_argument('-sys', help='system name - used to fetch gro and xtc filenames, and nprot, and nframes', choices=files_dict.keys(), default='none')
	parser_data.add_argument('-p', help='protein name - used to fetch nresm and nmono and protseq', choices=protein_info.keys(), default='Kir2.2_3spi')
	parser_data.add_argument('-reslist',default='all', help='''reslist should be input in the format:  1:4,6:9,300,310:333  \n(no spaces, residue ranges are inclusive,  \n
	Residue numbering starts at 0 - eg. residue 1 in your gro file is residue 0 when inputting here.)\n
	OR ALL residues can be specified using: all''')
	parser_data.add_argument('-syspath', help='system pathname - used to fetch bs info files - can be prespecified in files_dict', default='none')
	parser_data.add_argument('-f',help='input gro file - can be pre-specified in files_dict', default = 'none')
	parser_data.add_argument('-x', help='List of xtc files - can be pre-specified in files_dict', nargs='+', default = 'none')
	parser_data.add_argument('-np',default=1, help='Number of proteins in system - can be pre-specified in files_dict')
	parser_data.add_argument('-lipid', help='Name of lipid', required=True)
	parser_data.add_argument('-nframes', type=int, default = 0, help='Number of frames in trajectory - can be pre-specified in files_dict')
	parser_data.add_argument('-stride', type=int, default = 1, help='Frame number intervals that xtc will be read at')
	parser_data.add_argument('-lp', default='headgroup', help='Part of lipid to consider in interactions: headgroup or phosphate', choices=lipid_particles.keys())
	parser_data.add_argument('-d', type=float, default=6.5, help='Interactions cutoff distance (Angstroms)')
	#parser.add_argument('-dt', type=int, default=1, help='Number of nanoseconds per frame (once "stride" has been taken into account) - can be pre-specified in files_dict')
#	parser.add_argument('-permonomer', action='store_true', help='process results separately for each monomer - for this analysis you would not normally specify this')
	parser_data.add_argument('-nresm', type=int, default = 0, help='number of residues in each prot monomer - can be pre-specified in protein_info')
	parser_data.add_argument('-nmono', type=int, default = 0, help='Number of monomers per protein - can be pre-specified in protein_info')
	parser_data.add_argument('-seq', type=str, default = '', help='Protein sequence - can be pre-specified in protein_info')
	parser_data.add_argument('-fl', help='label to add to filenames - will automatically include sys and lipid', default = '')
	parser_data.add_argument('-fitting', type=str, default = 'double', help='type of fitting used to get bs info')
	parser_data.add_argument('-annotations', default = [], nargs='+', help='give list of bs annotations to focus on - if not specified will take from bs_info file')
	parser_data.set_defaults(func=process_data)

	parser_plot = subparsers.add_parser('plot', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_plot.add_argument('-sys', help='system name - used to fetch gro and xtc filenames, and nprot, and nframes', choices=files_dict.keys(), default='none')
	parser_plot.add_argument('-p', help='protein name - used to fetch nresm and nmono and protseq', choices=protein_info.keys(), default='Kir2.2_3spi')
	parser_plot.add_argument('-data', help='filename to read data from - if not specified will set to n_contacts_<file-label>.txt and if the -fl flag is not set, will set to n_contacts_<sys>_<lipid>.txt ', default = '')
	parser_plot.add_argument('-syspath', help='system pathname - used to fetch bs info files - can be prespecified in files_dict', default='none')
	parser_plot.add_argument('-lipid', help='Name of lipid', required=True)
	parser_plot.add_argument('-nresm', type=int, default = 0, help='number of residues in each prot monomer - can be pre-specified in protein_info')
	parser_plot.add_argument('-nmono', type=int, default = 0, help='Number of monomers per protein - can be pre-specified in protein_info')
	parser_plot.add_argument('-fl', help='label to add to filenames - will automatically include sys and lipid', default = '')
	parser_plot.add_argument('-fitting', type=str, default = 'double', help='type of fitting used to get bs info')
	parser_plot.add_argument('-annotations', default = [], nargs='+', help='give list of bs annotations to focus on - if not specified will take from bs_info file')
	parser_plot.set_defaults(func=plot_data)

	args = parser.parse_args()
	args.func(args)
