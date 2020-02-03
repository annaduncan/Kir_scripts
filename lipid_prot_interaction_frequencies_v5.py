'''Scheme of the script:
1. identify residue IDs for all repeats of protein (unecessary if only 1 protein obv)
?1a.  Create trajectory consisting of only relvant beads?
2. For each protein residue (or bead of residue) of interest:
    a. Calculate distances between lipid beads and protein residue
    b. If distance less than 0.7 and  add to total for the protein
    
Changes to v5 from v4 - re-write structure of count_frequencies algorithm
'''
import sys
import re
import numpy
import time
import MDAnalysis
from MDAnalysis.core.parallel.distances import distance_array
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cmx
import matplotlib.colors as colors

import ast
import argparse

from interaction_freq_core import protein_info, get_real_resnums_strlist, bondini_residuesymbol_size, martini_residuesymbol_size

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

def get_protresname_array():
	pass # could add this function later, to make barplots functions a bit more generic/smoother

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
#edge_colors_dict = { 'PM': 'orange','P2subPE': 'blue','darkgrey'#edge_colors norms PM, orange
#edge_colors norms PM, orange
#python /sansom/s105/bioc1280/Simulations/Tools/lipid_prot_interactionfreq/lipid_prot_interaction_frequencies_v5.py plot_ff -reslist all -i /sansom/n42/bioc1280/Kir22_3SPI_PM_3x3_ARCHER/prot_lipid_interactions/frequencies_dictionary_residuesall_lipidPOPSheadgroup_v2_5 /sansom/n42/bioc1280/Kir22_3SPIsys3_noGm3_3x3/prot_lipid_interactions/frequencies_dictionary_residuesall_lipidPOPSheadgroup_v2_5 /sansom/n42/bioc1280/Kir22_3SPIsys3_noCHOL3_3x3/prot_lipid_interactions/frequencies_dictionary_residuesall_lipidPOPSheadgroup_v2_5 /sansom/n42/bioc1280/Kir2_2_3SPI_PIP2subPE_3x3/prot_lipid_interactions/frequencies_dictionary_residuesall_lipidPOPSheadgroup_v2_5 /sansom/n42/bioc1280/Kir22_3SPI_noP2noGM3_3x3_ARCHER/prot_lipid_interactions/frequencies_dictionary_residuesall_lipidPOPSheadgroup_v2_5 /sansom/n42/bioc1280/Kir22_3SPIsys3_noP2noGm3noCHOL2_3x3/prot_lipid_interactions/frequencies_dictionary_residuesall_lipidPOPSheadgroup_v2_5 -p Kir2.2_3spi -lipid POPS POPS POPS POPS POPS POPS -s PM_3x3 noGM3_3x3 noCHOL_3x3 noP2_3x3 noGM3P2_3x3 noP2GM3CHOL_3x3

### main function to count interactions and produce interactions file ###
def count_frequencies(gro_file, trr_file, prot_seq, protein_residue_list, nrepeats, lipid, stride=1, lipid_part='headgroup', protein_centre='centroid', protein_centre_cutoff= 60, cutoff=6.5):

	universe = MDAnalysis.Universe(gro_file, trr_file)

	protein_res_total = len(prot_seq)
	protein_residue_dictionary = find_prot_residues(protein_res_total, protein_residue_list, nrepeats)
	prot_split = make_split_list_single_prot(prot_seq)
	#print prot_split
	
	lipid_selection = 'resname {} and ('.format(lipid)
	for bead in lipid_particles[lipid_part][lipid]:
		lipid_selection += 'name {} or '.format(bead)
	lipid_selection = lipid_selection[:-4] + ')'
	
	lipid_rep_selection = 'resname {} and name {}'.format(lipid, lipid_particles[lipid_part][lipid][0]) # ie. just choose one bead

	lipids = universe.selectAtoms(lipid_selection)
	n_lipid_beads = len(lipid_particles[lipid_part][lipid])
	n_lipids = lipids.numberOfAtoms() / n_lipid_beads
	lipid_reps = universe.selectAtoms(lipid_rep_selection)

    #initialise protein-lipid interactions frequency list
	proteinres_lipid_interactions = numpy.array([0 for i in protein_residue_list])

	startTime = time.time()    
	print 'Here we go...'
	frame = 0
	for ts in universe.trajectory[::stride]:
		for i in range(nrepeats):
			single_prot = universe.segments[0][range(i*protein_res_total,(i+1)*protein_res_total)]
			# find protein centroid - or pick out residue to represent protein position
			if protein_centre == 'centroid':
				single_prot_cent = numpy.array([single_prot.centroid()])
			elif (protein_centre) == int:
				#pick out BB of specified residue
				single_prot_cent = universe.segments[0][i*protein_res_total + protein_centre-1][0] # -1 is because res numbers are zero-indexed
			else:
				print 'Error: protein_centre should either be an integer residue number, or "centroid"'
				return None
			# find lipids within 60 A of prot centroid and pick those out from lipids selection
			close = distance_array(single_prot_cent, lipid_reps.coordinates(), ts.dimensions) < protein_centre_cutoff
			lipids_close_indices = []
			for index in numpy.nonzero(close)[1]: # numpy.nonzero gives a tuple of arrays - the second gives the indices of lipid residues that are 'close' to the protein
				lipids_close_indices += range(index*n_lipid_beads, (index+1)*n_lipid_beads) # convert residue IDs to atoms IDs for 'close' lipids
			lipids_close = lipids[lipids_close_indices]
			n_lipids_close = sum(close.flatten())
			# now look at prot-lipid interaction on per residue level
			if n_lipids_close == 0:
				continue # ie. nothing is added to proteinres_lipid_interactions
			else:
				all_dists = distance_array(single_prot.coordinates(), lipids_close.coordinates(), ts.dimensions)
				protein_lipid_dist_perresidue_all = numpy.array([[x.min() for x in numpy.split(lip, prot_split, axis=0)] for lip in numpy.split(all_dists, n_lipids_close, axis=1)])
				#print protein_lipid_dist_perresidue_all[:,17]
				interactions = protein_lipid_dist_perresidue_all <= cutoff
				proteinres_lipid_interactions += numpy.sum(interactions, axis = 0) # sum over all lipids for each res
			#update = '\rProtein {}/{} took {:3f} s'.format(i, nrepeats, time.time()-startTime)
			#print update,
			#sys.stdout.flush()
			#sys.stdout.write(update)
			#startTime = time.time()
		frame += 1
		print 'Frame {} took {:3f} s'.format(frame, time.time()-startTime)
		startTime = time.time()
	proteinres_lipid_interactions_dict = dict( zip(protein_residue_list, proteinres_lipid_interactions) )
	print proteinres_lipid_interactions_dict
	return proteinres_lipid_interactions_dict

### functions to produce figures from data ###
def plot_frequencies(protein_residue_list, proteinres_lipid_interactions):
    plt.bar(protein_residue_list, [proteinres_lipid_interactions[res] for res in protein_residue_list])
    plt.show()
            
def plot_frequencies_fancy(protein_residue_list, proteinres_lipid_interactions_dict, prot, sys_list, cutoff=True, cutoff_value=0.025):
	protein_residue_array = numpy.array(protein_residue_list)
	actual_res_numbers =  protein_residue_array + 40
	actual_res_namenums = numpy.core.defchararray.add(numpy.array(list(protein_info[prot]['protein_seq'])), get_real_resnums_strlist(prot))
	def percent_interactions(proteinres_lipid_interactions):
		total_interactions = sum(proteinres_lipid_interactions.values())
		interactions_list = []
		for res in protein_residue_list:
			interactions_list.append(float(proteinres_lipid_interactions[res]))
		interactions_array = numpy.array(interactions_list)
		percent_interactions = interactions_array*100/total_interactions
		return percent_interactions
	percent_interactions_dict = {}
	lipid_list = sorted(proteinres_lipid_interactions_dict.keys())
	for lipid in lipid_list:
		percent_interactions_dict[lipid] = percent_interactions(proteinres_lipid_interactions_dict[lipid])
	if cutoff:
		# graph will include ANY residue that has more than <<cutoff>> interactions for ANY lipid
		over_cutoff = numpy.zeros(shape=len(protein_residue_list), dtype=bool)
		for lipid in lipid_list:
			over_cutoff += percent_interactions_dict[lipid] >= cutoff_value*100
		for lipid in lipid_list:
			percent_interactions_dict[lipid] = percent_interactions_dict[lipid][over_cutoff]
		actual_res_namenums = actual_res_namenums[over_cutoff]
	fig, ax = plt.subplots()
	#rects=ax.barh(range(len(percent_interactions)), percent_interactions, align='center' , color='yellow')
	#ax.invert_xaxis()
	ax.invert_yaxis()
	width = 0.7/len(lipid_list)
	#edge_colors = ['orange','blue','darkgrey']
	edge_colours = [sys_colors[sys] for sys in sys_list]
	print edge_colours
	i=0
	lipid_rects = {}
	for lipid in lipid_list:
		if re.match(r'PI.*', lipid):
			#barcolor = 'yellow'
			#maxpercent = max(percent_interactions_dict[lipid])
			lipid_cmap = plt.get_cmap('hot_r')
		elif re.match(r'CHOL.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('cyanmap', ['white', 'darkcyan'])
		elif re.match(r'POPS.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('bluemap', ['white','darkblue'])
		elif re.match(r'GM3.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('magentamap', ['white', 'darkmagenta'])
		elif re.match(r'POPC.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('greymap', ['white','black'])
		elif re.match(r'POPE.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('greymap', ['white','darkgrey'])
		elif re.match(r'PPCS.*', lipid):
			lipid_cmap=colors.LinearSegmentedColormap.from_list('greenmap', ['white','green'])
		maxpercent = max(max(percent_interactions_dict[lipid]) for lipid in lipid_list)
		cNorm = colors.Normalize(vmin=0,vmax=maxpercent)
		scalar_map = cmx.ScalarMappable(norm = cNorm, cmap = lipid_cmap)
		barcolor=scalar_map.to_rgba(percent_interactions_dict[lipid])
		#lipid_rects[lipid] = ax.barh(numpy.arange(len(percent_interactions_dict[lipid]))+width*i, percent_interactions_dict[lipid], width, align='center' , color=barcolor, ec=edge_colors[i])
		lipid_rects[lipid] = ax.barh(numpy.arange(len(percent_interactions_dict[lipid]))+width*i, percent_interactions_dict[lipid], width, align='center' , color=edge_colours[i])#, ec=barcolor)
		i+=1
	#if len(proteinres_lipid_interactions_dict.keys()) > 1:
	#	ax.legend( lipid_rects.values(), lipid_rects.keys() )
	ax.set_ylabel('Residue')
	ax.set_xlabel('% of total {} contacts'.format(lipid_list[0]))
	ax.set_yticks( range(len(actual_res_namenums)))
	ax.set_yticklabels(actual_res_namenums, rotation='horizontal')
	plt.axis('tight')
	fig.set_size_inches(3.7,8.3)
	plt.savefig('prot_{}_int_v6.svg'.format(''.join(lipid_list)), orientation='portrait', bbox_inches='tight')
	plt.savefig('prot_{}_int_v6.png'.format(''.join(lipid_list)), orientation='portrait', bbox_inches='tight')
	#plt.show()

def show_frequencies_on_structure(nmonomers, interactions, gro_file, outfile, sensitive=False):
	#normalise interactions - not really necessary - can do this w/ VMD colorscales
	u = MDAnalysis.Universe(gro_file)
	protein=u.selectAtoms("protein")
	nres = len(protein.residues)
	nresmono = nres/nmonomers
	total_interactions = sum(interactions.values())
	protein_residue_list = sorted(interactions.keys())
	interactions_array = numpy.array([float(interactions[res]) for res in protein_residue_list])
	if sensitive:
		percent_interactions = interactions_array*10000/total_interactions # to make smaller numbers visible (for when total_interactions is v large)
	else:
		percent_interactions = interactions_array*100/total_interactions
	for i in range(nresmono):
		for j in range(nmonomers):
			protein.residues[j*nresmono + i].set_bfactor(percent_interactions[i])
	with MDAnalysis.Writer(outfile, numatoms=u.atoms.numberOfAtoms()) as PDB:
		PDB.write(u.atoms)

# functions for reading input - residue list and interactions file
def make_res_list(res_string, protname):
	# residue list should be input in the format 1:4,6:9,300,310:333
	# the following decodes that format into a straight list
	if res_string == 'all':
		res_list = numpy.arange(len(list(protein_info[protname]['protein_seq']))) + 1
	else:
		res_list_ranges = [x.split(':') for x in res_string.split(',')]
		res_list = []
		for res_range in res_list_ranges:
			res_list += range(int(res_range[0]), int(res_range[-1])+1 ) 
	#print 'res_list: ', res_list
	return res_list
def make_interactions_dict(dictionary_filename):
	f = open(dictionary_filename)
	dictionary_string = f.readlines()[2].rstrip()
	proteinres_lipid_interactions = ast.literal_eval(dictionary_string)
	return proteinres_lipid_interactions

### Function to sum data for all (or a selection of) lipids:
def sum_lipid_ints(interactions_dict_list):
	total_interactions = {}
	for interactions_dict in interactions_dict_list:
		for res in interactions_dict.keys():
			if res not in total_interactions:
				total_interactions[res] = interactions_dict[res]
			else:
				total_interactions[res] += interactions_dict[res]
	return total_interactions

if __name__ == "__main__":
	# sub-command funcs for the parser - ie. workflows available when running script:
	def get_ints(args):
		res_list = make_res_list(args.reslist, args.p)
		protein_seq = list(protein_info[args.p]['protein_seq'])
		if args.c != 'centroid':
			args.c = int(args.c) # this should return an error if the argument is not integer-like
			if args.c < len(protein_seq):
				#check that the residue number is valid
				print 'Error: residue number selected is larger than the number of residues in the protein sequence'
				return None
		proteinres_lipid_interactions = count_frequencies(args.f, args.xtc, protein_seq, res_list, args.nrepeats, args.lipid, args.stride, args.lipid_part, args.c, args.cd, args.cutoff)
		f=open('frequencies_dictionary_residues{}_lipid{}{}_v2_5'.format(args.reslist,args.lipid,args.lipid_part), 'w')
		f.write('Generated using the command:\npython {}\n'.format(' '.join(sys.argv)))
		f.write(str(proteinres_lipid_interactions)+'\n')
		f.close()

	def pff(args):
		res_list = make_res_list(args.reslist, args.p)
		proteinres_lipid_interactions_dict = {}
		for i in range(len(args.lipid)):
			if args.lipid[i] not in proteinres_lipid_interactions_dict:
				proteinres_lipid_interactions_dict[args.lipid[i]] = make_interactions_dict(args.i[i])
				print '{} added with lipid name {}'.format(args.i[i],args.lipid[i])
			elif args.lipid[i]+'x' not in proteinres_lipid_interactions_dict:
				proteinres_lipid_interactions_dict[args.lipid[i]+'x'] = make_interactions_dict(args.i[i])
				print '{} added with lipid name {}'.format(args.i[i],args.lipid[i]+'x')
			elif args.lipid[i]+'xx' not in proteinres_lipid_interactions_dict:
				proteinres_lipid_interactions_dict[args.lipid[i]+'xx'] = make_interactions_dict(args.i[i])
				print '{} added with lipid name {}'.format(args.i[i],args.lipid[i]+'xx')
			elif args.lipid[i]+'xxx' not in proteinres_lipid_interactions_dict:
				proteinres_lipid_interactions_dict[args.lipid[i]+'xxx'] = make_interactions_dict(args.i[i])
				print '{} added with lipid name {}'.format(args.i[i],args.lipid[i]+'xxx')
			elif args.lipid[i]+'xxxx' not in proteinres_lipid_interactions_dict: 
				proteinres_lipid_interactions_dict[args.lipid[i]+'xxxx'] = make_interactions_dict(args.i[i])
				print '{} added with lipid name {}'.format(args.i[i],args.lipid[i]+'xxxx')
			else: # this will do for 6 systems but need to update if comparing more...
				proteinres_lipid_interactions_dict[args.lipid[i]+'xxxxx'] = make_interactions_dict(args.i[i])
				print '{} added with lipid name {}'.format(args.i[i],args.lipid[i]+'xxxxx')
				
		#print proteinres_lipid_interactions_dict
		proteinres_lipid_interactions = plot_frequencies_fancy(res_list, proteinres_lipid_interactions_dict, args.p, args.s)

	def structure(args):
		interactions = make_interactions_dict(args.i)
		show_frequencies_on_structure(args.m, interactions, args.f, args.o, args.s)
	
	def sum_ints(args):
		interactions_dict_list = []
		for i in range(len(args.i)):
			interactions_dict_list.append(make_interactions_dict(args.i[i]))
		total_interactions = sum_lipid_ints(interactions_dict_list)
		f=open('frequencies_dictionary_sumlipids_v2_5'.format(), 'w')
		f.write('Generated using the command:\npython {}\n'.format(' '.join(sys.argv)))
		f.write(str(total_interactions)+'\n')
		f.close()
		
		
	parser = argparse.ArgumentParser(description='Analyse prot-lipid interactions at the residue level', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = parser.add_subparsers()

	parser_ints = subparsers.add_parser('get_ints', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_ints.add_argument('-f',help='input gro file', required=True)
	parser_ints.add_argument('-x', '--xtc', help='List of xtc files', nargs='+', required=True)
	parser_ints.add_argument('-reslist',default='all', help='''reslist should be input in the format:  1:4,6:9,300,310:333  (no spaces, residue ranges are inclusive) 
	Residue numbering starts at 1)
	OR ALL residues can be specified using: all''')
	parser_ints.add_argument('-p', help='Protein name', choices=protein_info.keys(), required=True)
	parser_ints.add_argument('-r', '--nrepeats', type=int, help='Number of protein MONOMERS in system', required=True)
	parser_ints.add_argument('-lipid', help='Name of lipid', choices=lipid_particles['headgroup'], required=True)
	parser_ints.add_argument('-stride', type=int, default = 1, help='Frame number intervals that xtc will be read at')
	parser_ints.add_argument('-lp', '--lipid_part', default='headgroup', help='Part of lipid to consider in interactions', choices=lipid_particles.keys())
	parser_ints.add_argument('-c', default='centroid', help='Part of protein to measure lipid distances from, for first approximation - can be either "centroid" for protein centroid, or an integer residue number (residue numbering starts from 1)')
	parser_ints.add_argument('-cd', type = float, default=80, help='Distance away from protein to select lipids that will undergo a closer inspection (a lower number speeds the script up, but too low may mean that some lipids aren\'t considered')
	parser_ints.add_argument('-d', '--cutoff', type=float, default=6.5, help='Interactions cutoff distance')
	parser_ints.set_defaults(func=get_ints)

	parser_pff = subparsers.add_parser('plot_ff', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_pff.add_argument('-reslist',default='all', help='''reslist should be input in the format:  1:4,6:9,300,310:333  (no spaces, residue ranges are inclusive - 
	Residue numbering starts at 1)
	OR ALL residues can be specified using: all''')
	parser_pff.add_argument('-i',help='Interaction file - can take multiple and display on one bar chart', nargs='+')
	parser_pff.add_argument('-s',help='list of systems (to give edge colour)', nargs='+')
	parser_pff.add_argument('-p', help='Protein name', choices=protein_info.keys(), required=True)
	parser_pff.add_argument('-lipid', default='PIP2', help='Name of lipid - list lipid multiple times if it is the same lipid but a different interaction file', choices=lipid_particles['headgroup'], nargs='+')
	parser_pff.set_defaults(func=pff)

	parser_struct = subparsers.add_parser('structure', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_struct.add_argument('-i', help='Interaction file')
	parser_struct.add_argument('-m',type=int, help='Number of monomers in one protein')
	parser_struct.add_argument('-f',help='input gro file')
	parser_struct.add_argument('-o',help='output gro filename')
	parser_struct.add_argument('-s',help='sensitive - amplifies all interactions in order to pick up the least frequent interactions', action='store_true', default=False)
	parser_struct.set_defaults(func=structure)
	
	parser_addints = subparsers.add_parser('add', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_addints.add_argument('-i', help='Interaction files', nargs = '+', required=True)
	parser_addints.set_defaults(func=sum_ints)
	
	args = parser.parse_args()
	args.func(args)

#	parser.add_argument('-i','--interaction_files', help='List of interaction files', nargs='+', type = lambda s: [item for item in s.split(' ')])
#	parser.add_argument('-s','--systems_list', help='List of system names given in interaction files list', nargs='+', type = lambda s: [item for item in s.split(' ')]
	

	
	
        
                  
        
         
         
