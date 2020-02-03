'''Scheme of the script:
1. identify residue IDs for all repeats of protein (unecessary if only 1 protein obv)
?1a.  Create trajectory consisting of only relvant beads?
2. For each protein residue (or bead of residue) of interest:
    a. Calculate distances between lipid beads and protein residue
    b. If distance less than 0.7 and  add to total for the protein
    
Changes to v5 from v4 - re-write structure of count_frequencies algorithm
13/March/2017 - change dt from int to float
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
import networkx as nx
import community

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

import time

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
	'headgroup' : {'PIP2': ['PO1','PO2','PO3','RP1','RP2','RP3'],
			'PI3': ['PO0','PO1','PO2','PO3','RP1','RP2','RP3'],
			'PSPI' : ['PO3','RP1','RP2','RP3'],
			'CHOL' : ['ROH', 'R1', 'R2', 'R3', 'R4', 'R5'], #'C1', 'C2'], # this is actually the whole CHOL},
			'GM3': ['B1A','B2A','B3A','INA','B1B','B2B','B3B','INB','B1C','B2C','B3C','INC','B4C','B5C'],
			'PPCS' : ['PO4','NC3'],
			'POPS' : ['PO4','CNO'],
			'POPC' : ['PO4','NC3'],
			'POPE' : ['PO4','NH3'],
			'MLCA' : ['PO41','PO42','GL0'],
			'MLCB' : ['PO41','PO42','GL0'],
			'CDL4' : ['PO41','PO42','GL0'],
			'CDL2' : ['PO41','PO42','GL0']},
	'phosphate' : {'PIP2': ['PO3','PO1','PO2'],
			'PI3': ['PO3','PO0','PO1','PO2'],
			'PSPI' : ['PO3'],
			'CHOL' : ['ROH'],
			'GM3': ['B1A','B2A','B3A','INA'], # This is obviously not a phosphate group, but it's in the position of one in the memb
			'PPCS' : ['PO4'],
			'POPS' : ['PO4'],
			'POPC' : ['PO4'],
			'POPE' : ['PO4'],
			'MLCA' : ['PO41','PO42'],
			'MLCB' : ['PO41','PO42'],
			'CDL4' : ['PO41','PO42'],
			'CDL2' : ['PO41','PO42']}
}

#protein_sites = {'Kir2.2_3spi':{'primary':numpy.array([[37, 39, 142, 145, 147, 148], [369, 371, 474, 477, 479, 480], [701, 703, 806, 809, 811, 812], [1033, 1035, 1138, 1141, 1143, 1144]]), 'secondary_K220':numpy.array([[843], [1175], [511], [179]]), 'secondary_K62':numpy.array([[21], [353], [685], [1017]])}} # these are the zero-indexed gro file residue numbers, not the actual res IDs

#protein_sites = {'Kir2.2_3spi':{'primary_R78':numpy.array([[37], [369], [701], [1033]]), 'primary_R80':numpy.array([[39], [371], [703], [1035]]), 'primary_K183':numpy.array([[142], [474], [806], [1138]]), 'primary_R186':numpy.array([[145], [477], [809], [1141]]), 'primary_K188':numpy.array([[147], [479], [811], [1143]]), 'primary_K189':numpy.array([[148], [480], [812], [1144]]), 'secondary_K220':numpy.array([[843], [1175], [511], [179]]), 'secondary_K62':numpy.array([[21], [353], [685], [1017]])}} # these are the zero-indexed gro file residue numbers, not the actual res IDs

#sites_index_colours = {'Kir2.2_3spi':{'primary':3, 'secondary_K220':2, 'secondary_K62':1}}
sites_index_colours = {'Kir2.2_3spi':{'primary_R78':3, 'primary_R80':4, 'primary_K183':5, 'primary_R186':6, 'primary_K188':7, 'primary_K189':8, 'secondary_K220':2, 'secondary_K62':1}}

def centroid_array_production_protein(protein_sel,num_protein_copies):
    dictionary_centroid_arrays = {}
    full_protein_coord_array = protein_sel.coordinates()
    list_individual_protein_coord_arrays = numpy.split(full_protein_coord_array,num_protein_copies)
    list_per_protein_centroids = [numpy.average(protein_coord_array,axis=0) for protein_coord_array in list_individual_protein_coord_arrays]
    dictionary_centroid_arrays['protein'] = numpy.array(list_per_protein_centroids)
    return dictionary_centroid_arrays


### main function to collate interaction distances and produce interaction distance covariances ###
def fetch_interactions(gro_file, trr_file, prot_seq, prot_name, nmonomers, lipid, prot_index_list, nframes, stride=1, lipid_part='headgroup', cutoff=65):
	universe = MDAnalysis.Universe(gro_file, trr_file)

	protein_res_total = len(prot_seq)
	#protein_residue_dictionary = find_prot_residues(protein_res_total, protein_residue_list, nrepeats=1)
	prot_split_monomer = make_split_list_single_prot(prot_seq)
	prot_split = []
	n_prot_atoms_permonomer = len(universe.segments[0][numpy.arange(protein_res_total)].atoms)
	print 'n_prot_atoms_permonomer', n_prot_atoms_permonomer
	for m in range(nmonomers):
		prot_split += list( numpy.array(prot_split_monomer) + m*(n_prot_atoms_permonomer) )
		if (m+1) < (nmonomers):
			prot_split += [(m+1)*(n_prot_atoms_permonomer)]

	
	lipid_selection = 'resname {} and ('.format(lipid)
	for bead in lipid_particles[lipid_part][lipid]:
		lipid_selection += 'name {} or '.format(bead)
	lipid_selection = lipid_selection[:-4] + ')'
	
	#lipid_rep_selection = 'resname {} and name {}'.format(lipid, lipid_particles[lipid_part][lipid][0]) # ie. just choose one bead

	lipids = universe.selectAtoms(lipid_selection)
	n_lipid_beads = len(lipid_particles[lipid_part][lipid])
	n_lipids = lipids.numberOfAtoms() / n_lipid_beads
	#lipid_reps = universe.selectAtoms(lipid_rep_selection)
	lipid_indices = lipids.residues.resids()
	
	startTime = time.time()    
	print 'Here we go...'
	frame = 0
	prot_lipid_dists = {}
	#sites_list = protein_sites[prot_name].keys()

	for protein_index in prot_index_list:
		prot_lipid_dists[protein_index] = {}
		prot_lipid_dists[protein_index] = numpy.zeros((protein_res_total*nmonomers, n_lipids, nframes),dtype=int)
	for ts in universe.trajectory[::stride]:
		if frame >= nframes:
			print 'Have reached maximum number of frames specified.  Stopping...'
			continue
		for protein_index in prot_index_list:
			## get min dist per lipid and residue
			single_prot = universe.segments[0][protein_res_total*protein_index*nmonomers + numpy.arange(protein_res_total * nmonomers) ]
			dists = distance_array(single_prot.coordinates(), lipids.coordinates(), ts.dimensions)
			min_dists_per_lipid = numpy.min(numpy.array(numpy.split(dists, n_lipids, axis=1)), axis = 2)
			#if frame == 0:
			#	print 'min_dists_per_lipid.shape', min_dists_per_lipid.shape
			split_per_res = numpy.split(min_dists_per_lipid, prot_split, axis=1)
			min_dists_perresidue = numpy.array([x.min(axis=1) for x in split_per_res] )
			#if frame == 0:
			#	print 'min_dists_perresidue.shape', min_dists_perresidue.shape
			prot_lipid_dists[protein_index][:,:,frame] = min_dists_perresidue
		if frame == 0:
			print 'Frame {} (fromtraj)or{} (hardcoded) took {:3f} s\r'.format(ts.frame, frame, time.time()-startTime)
		frame += 1
		startTime = time.time()
	return prot_lipid_dists, lipid_indices

def sum_over_monomers(interaction_matrix, nmonomers):
    summed_matrix = numpy.sum(numpy.array(numpy.split(numpy.sum(numpy.array(numpy.split(interaction_matrix,nmonomers,axis=1)), axis=0), nmonomers,axis=0)), axis=0)
    return summed_matrix
    
def interaction_covariances(prot_lipid_dists, prot_seq, nmonomers, n_lipids, nframes, measure = 'correlation', compare='cutoff', cutoff=7):
	#distance_covariances = numpy.zeros( len(prot_seq) * nmonomers )
	#flattened_prot_lipid_dists = prot_lipid_dists[0].reshape(len(prot_seq)*nmonomers, n_lipids*nframes)
	nprots = len(prot_lipid_dists)
	nres = len(prot_seq)
	print nprots, nres, nmonomers
	all_prots_distances = numpy.zeros( (nres*nmonomers, nprots*n_lipids*nframes) )
	for prot in range(nprots):
		all_prots_distances[:,prot*n_lipids*nframes:(prot+1)*n_lipids*nframes] = prot_lipid_dists[prot].reshape(nres*nmonomers, n_lipids*nframes)
	all_prots_count  = all_prots_distances < cutoff
	all_prots_count_allmonomers = numpy.sum(all_prots_count, axis=1)
	print 'all_prots_count_allmonomers.shape ', all_prots_count_allmonomers.shape
	#print 'all_prots_count_allmonomers ', all_prots_count_allmonomers
	single_prot_count = numpy.sum(numpy.array(numpy.split(all_prots_count_allmonomers,nmonomers)), axis=0)
	#print all_prots_count
	single_prot_count = single_prot_count *300. /numpy.max(single_prot_count)
	print 'all_prots_count_allmonomers.shape', all_prots_count_allmonomers.shape, 'single_prot_count.shape', single_prot_count.shape
	#print 'single_prot_count.max', numpy.max(single_prot_count)
	#print single_prot_count
	#near_residues = all_prots_distances[all_prots_distances < 7]
	near_residues_all_frames = numpy.any(all_prots_distances < cutoff, axis = 1)
	near_residues_all_frames_avg_monomer = numpy.sum(numpy.array(numpy.split(near_residues_all_frames, nmonomers)),axis=0)
	if compare == 'cutoff':
		all_near_residues_all_frames = (all_prots_distances < cutoff)*1
		#print all_near_residues_all_frames
		print 'all_near_residues_all_frames.shape', all_near_residues_all_frames.shape
		#print all_near_residues_all_frames[near_residues_all_frames]
		print 'type(all_near_residues_all_frames)', type(all_near_residues_all_frames)
		if measure == 'covariance':
			near_distance_covariances = numpy.cov(all_near_residues_all_frames[near_residues_all_frames])
		elif measure == 'correlation':
			near_distance_covariances = numpy.corrcoef(all_near_residues_all_frames[near_residues_all_frames])
#		distance_covariances = near_distance_covariances
	print 'near_residues_all_frames.shape', near_residues_all_frames.shape
	print 'all_prots_distances.shape should be nprots*nres*nmonomers by nlipids*nframes: ', all_prots_distances.shape
	if compare == 'dists':
		if measure == 'covariance':
			near_distance_covariances = numpy.cov(all_prots_distances[near_residues_all_frames])
		elif measure == 'correlation':
			near_distance_covariances = numpy.corrcoef(all_prots_distances[near_residues_all_frames])
	n_near_residues = near_distance_covariances.shape[0]
	print 'near_distance_covariances.shape', near_distance_covariances.shape
	distance_covariances = numpy.zeros((nres*nmonomers, nres*nmonomers)) # nprots*n_lipids*nframes))
	indexes_to_put_into_whole_matrix = numpy.tile(near_residues_all_frames, (nres*nmonomers,1)) * numpy.tile(near_residues_all_frames, (nres*nmonomers,1)).T
	print 'near_residues_all_frames.shape', near_residues_all_frames.shape
	print 'indexes_to_put_into_whole_matrix.shape', indexes_to_put_into_whole_matrix.shape
	distance_covariances[indexes_to_put_into_whole_matrix] = near_distance_covariances.flatten()
	print 'Finished calculating covariances'
	distance_covariances = numpy.sum(numpy.array(numpy.split(numpy.sum(numpy.array(numpy.split(distance_covariances,nmonomers,axis=1)), axis=0), nmonomers,axis=0)), axis=0)
	print 'distance_covariances.shape', distance_covariances.shape
	#near_distance_covariances = distance_covariances[near_residues_all_frames_avg_monomer][:,near_residues_all_frames_avg_monomer]
	#near_distance_covariances = numpy.sum(numpy.array(numpy.split(numpy.sum(numpy.array(numpy.split(near_distance_covariances,4,axis=1)), axis=0), 4,axis=0)), axis=0)
	return distance_covariances, near_distance_covariances, near_residues_all_frames, single_prot_count

def graph_graph(graph_object,outputfilename, interaction_strength=numpy.array([0]), layout='spring', show_labels=True, node_colour='r', first_res_num=41):
	fig, ax = plt.subplots(nrows = 1)
	if layout == 'spring':
		pos=nx.spring_layout(graph_object)
		print pos
	elif layout == 'circular':
		pos=nx.circular_layout(graph_object)
		print pos
	if interaction_strength.shape == (1,):
		nx.draw_networkx_nodes(graph_object,pos)
	else:
		# make sure interaction strength is in the same order as graph.nodes()
		graph_nodes = graph_object.nodes()
		ordered_interaction_strength = [interaction_strength[sorted(graph_nodes).index(node)] for node in graph_nodes]
		nx.draw_networkx_nodes(graph_object,pos, node_size = ordered_interaction_strength)
		print graph_object.nodes(), interaction_strength
	if show_labels:
		lable_dict = dict(zip(graph_object.nodes(), numpy.array(graph_object.nodes())+first_res_num) )
#		print lable_dict
		nx.draw_networkx_labels(graph_object,pos,font_size=8, labels = lable_dict, font_color='gray' )
		print lable_dict
	weights = [data[2]['weight'] for data in graph_object.edges(data=True)]
	nx.draw_networkx_edges(graph_object,pos, width=weights)
	#nx.draw_networkx_edges(connected_component,pos)
	plt.axis('off')
	plt.savefig(outputfilename)
	plt.savefig(outputfilename[:-4]+'.pdf')
	plt.close()

def write_node_file(node_list, int_strength, filename, edge_weights=''):
	f = open(filename, 'w')
	f.write('# node list, followed by interaction strength for each node\n')
	f.write(str(node_list)+'\n')
	f.write(str(int_strength)+'\n')
	f.write(str(edge_weights)+'\n')
	f.close()

def network_based_on_correlation_analysis(distance_covariances, residue_interaction_strength_count, cov_cutoff=1.0, residue_interaction_network_connected_component_refined=False, residue_interaction_network_connected_component_unweighted=False, residue_interaction_network_clique_unweighted=False, unweighted_modularity_network=False, first_resnum=41):
	residue_network_raw = nx.Graph(distance_covariances)
	graph_graph(residue_network_raw, 'residue_interaction_network.svg',layout='circular', show_labels=False, interaction_strength=residue_interaction_strength_count, first_res_num=first_resnum)

	#residue_network_refined = nx.Graph(distance_covariances>1.])
	connected_graphs = list(nx.connected_component_subgraphs(residue_network_raw))
	for i,connected_component in enumerate(connected_graphs):
		if len(connected_component.nodes()) == 1:
			continue
		int_strength = residue_interaction_strength_count[connected_component.nodes()]
		graph_graph(connected_component,'residue_interaction_network_connected_component_{}.svg'.format(i), interaction_strength=int_strength, first_res_num=first_resnum)
		#graph_graph(connected_component,'residue_interaction_network_connected_component_{}.svg'.format(i), first_resnum=first_resnum)

	if residue_interaction_network_connected_component_refined:
		refined_network =numpy.copy(distance_covariances)
		refined_network[distance_covariances < cov_cutoff] = 0
		residue_network_refined = nx.Graph(refined_network)
		connected_refined_graphs = list(nx.connected_component_subgraphs(residue_network_refined))
		#print 'connected_refined',[graph.nodes() for graph in connected_refined_graphs]
		for i,connected_component in enumerate(connected_refined_graphs):
			if len(connected_component.nodes()) == 1:
				continue
			graph_graph(connected_component,'residue_interaction_network_connected_component_refined_{}.svg'.format(i), interaction_strength=int_strength, first_res_num=first_resnum)
	
	if residue_interaction_network_connected_component_unweighted:
		connected_refined_graphs = list(nx.connected_component_subgraphs(residue_network_refined))
		#print 'connected_refined',[graph.nodes() for graph in connected_refined_graphs]
		for i,connected_component in enumerate(connected_refined_graphs):
			if len(connected_component.nodes()) == 1:
				continue
			int_strength = residue_interaction_strength_count[connected_component.nodes()]
			graph_graph(connected_component,'residue_interaction_network_connected_component_unweighted_{}.svg'.format(i), interaction_strength=int_strength, first_res_num=first_resnum)

	if residue_interaction_network_clique_unweighted:
		unweighted_graph_array = numpy.copy(distance_covariances)
		unweighted_graph_array[distance_covariances < cov_cutoff] = 0
		unweighted_graph_array[distance_covariances >= cov_cutoff] = 1
		unweighted_graph = nx.Graph(unweighted_graph_array)
		cliques = list(nx.find_cliques(unweighted_graph))
		#print cliques
		for i,clique in enumerate(cliques):
			if len(clique) == 1:
				continue
			clique_graph = nx.subgraph(residue_network_raw, clique)
			int_strength = residue_interaction_strength_count[clique]
			graph_graph(clique_graph,'residue_interaction_network_clique_unweighted_{}.svg'.format(i), layout='circular', interaction_strength=int_strength, first_res_num=first_resnum)

	
	part = community.best_partition(residue_network_raw, weight='weight')
	mod = community.modularity(part,residue_network_raw)
	values = [part.get(node) for node in residue_network_raw.nodes()]
	for value in range(max(values)):
		node_list = [k for k,v in part.items() if v == value]
		if len(node_list) == 1:
			continue
		int_strength = residue_interaction_strength_count[node_list]
		subcommunity = nx.subgraph(residue_network_raw, node_list)
		graph_graph(subcommunity,'modularity_network_{}.svg'.format(value), interaction_strength=int_strength, first_res_num=first_resnum)
		write_node_file(node_list, int_strength, 'modularity_network_{}.txt'.format(value),distance_covariances[node_list,:][:,node_list])

	if unweighted_modularity_network:
		unweighted_graph_array = numpy.copy(distance_covariances)
		unweighted_graph_array[distance_covariances < 0.01] = 0
		unweighted_graph_array[distance_covariances >= 0.01] = 1
		unweighted_graph = nx.Graph(unweighted_graph_array)
		part_unweighted = community.best_partition(unweighted_graph)
		mod = community.modularity(part_unweighted,unweighted_graph)
		values = [part_unweighted.get(node) for node in unweighted_graph.nodes()]
		for value in range(max(values)):
			node_list = [k for k,v in part_unweighted.items() if v == value]
			if len(node_list) == 1:
				continue
			int_strength = residue_interaction_strength_count[node_list]
			subcommunity = nx.subgraph(residue_network_raw, node_list)
			graph_graph(subcommunity,'unweighted_modularity_network_{}.svg'.format(value), interaction_strength=int_strength, first_res_num=first_resnum)

def covariance_or_correlation_matrix_image(distance_covariances, measure = 'correlation', compare='cutoff'):
	fig, ax = plt.subplots(nrows = 1)
	cb = ax.imshow(distance_covariances, cmap='BuPu', interpolation='none')
	fig.colorbar(cb)
	if compare == 'cutoff':
		if measure == 'covariance':
			plt.savefig('cov_cutoff.svg')
		elif measure == 'correlation':
			plt.savefig('corrcoef_cutoff.svg')
	elif compare == 'dists':
		if measure == 'covariance':
			plt.savefig('cov.svg')
		elif measure == 'correlation':
			plt.savefig('corrcoef.svg')
	plt.close()

def covariance_or_correlation_matrix_image_closeonly(distance_covariances, near_residues_all_frames, Nterm_start, measure = 'correlation', compare='cutoff'):
	fig, ax = plt.subplots(nrows = 1)
	cb = ax.imshow(distance_covariances, cmap='BuPu', interpolation='none')
	labels = numpy.nonzero(near_residues_all_frames)[0] + Nterm_start
	plt.xticks(range(len(labels)), labels, fontsize=8)
	fig.colorbar(cb)
	if compare == 'cutoff':
		if measure == 'covariance':
			plt.savefig('cov_cutoff_close.svg')
		elif measure == 'correlation':
			plt.savefig('corrcoef_cutoff_close.svg')
	elif compare == 'dists':
		if measure == 'covariance':
			plt.savefig('cov_close.svg')
		elif measure == 'correlation':
			plt.savefig('corrcoef_close.svg')
	plt.close()

def dendrogram_image(distance_covariances, near_residues_all_frames, Nterm_start, compare='cutoff'):
	z = linkage(distance_covariances, method='centroid')
	z_labels = numpy.nonzero(near_residues_all_frames)[0] + Nterm_start
	fig = plt.figure()
	dn = dendrogram(z, 20, 'level', labels=z_labels)
	if compare == 'dists':
		plt.savefig('dendrogram.svg')
	elif compare == 'cutoff':
		plt.savefig('dendrogram_cutoff.svg')
	plt.close()
		
def clusters_on_structure(structure, distance_covariances, near_residues_all_frames, Nterm_start, compare='cutoff'):
	z = linkage(distance_covariances, method='centroid')
	clusters = fcluster(z, 50, criterion='distance')
	print clusters
	u = MDAnalysis.Universe(structure)
	protein=u.selectAtoms("protein")
	nres = len(protein.residues)
	for i in range(nres):
		protein.residues[i].set_bfactor(clusters[i])
	with MDAnalysis.Writer('clusters_on_struct.pdb', numatoms=u.atoms.numberOfAtoms()) as PDB:
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

def make_prot_index_list(prot_index_string):
	prot_ranges = [x.split(':') for x in prot_index_string.split(',')]
	prot_index_list = []
	for prot_range in prot_ranges:
		prot_index_list += range(int(prot_range[0]), int(prot_range[-1])+1 )
	return prot_index_list

def make_interactions_dict(dictionary_filename):
	f = open(dictionary_filename)
	dictionary_string = f.readlines()[2].rstrip()
	proteinres_lipid_interactions = ast.literal_eval(dictionary_string)
	return proteinres_lipid_interactions

#def distances_image_single_graph(prot_lipid_dists_image_dict, nframes, lipid, nmonomers, lipid_indices, sites_list, sites_index_colours_dict, output_filename, dt, show_all, monomers_sep):
def distances_image_single_graph(prot_lipid_dists_image_dict, nframes, lipid_list, nmonomers, lipid_indices_dict, sites_list, sites_index_colours_dict, output_filename, dt, show_all, monomers_sep):
	if monomers_sep:
		nmonomers_todisplay = 1
	else:
		nmonomers_todisplay = nmonomers
	lipids_to_show = {}
	prot_dists_image = {}
	nsites = len(sites_list)
	n_lipid_types = len(lipid_list)
	for lipid in lipid_list:
		if show_all:
			lipids_to_show[lipid] = numpy.array([True]*(prot_lipid_dists_image_dict[lipid][sites_list[0]].shape[0]))
		else:
			#lipids_or_array = prot_lipid_dists_image_dict[sites_list[0]].any(axis=1)
			#for site in sites_list[1:]:
			#	lipids_or_array += prot_lipid_dists_image_dict[site].any(axis=1)
			# show any lipid that has an interaction
			#lipids_to_show = lipids_or_array > 0
			lipids_or_array = numpy.sum(prot_lipid_dists_image_dict[lipid][sites_list[0]], axis=1)
			for site in sites_list[1:]:
				lipids_or_array += numpy.sum(prot_lipid_dists_image_dict[lipid][site], axis=1)
			#print 'lipids_or_array.shape', lipids_or_array.shape
			lipids_to_show[lipid] = lipids_or_array>10
		n_lipids = numpy.sum(lipids_to_show[lipid])
		prot_dists_image[lipid] = numpy.zeros((nsites*n_lipids, nframes), dtype=int)
		for site in sites_list:
			site_index = sites_index_colours_dict[site]
			prot_dists_image[lipid][[nsites*i+(site_index - 1) for i in range(n_lipids)],] = (prot_lipid_dists_image_dict[lipid][site]*site_index)[lipids_to_show[lipid]]
	if nframes < 100:
		figsize_x = 8
	else:
		figsize_x = 8*nframes/400.
	fig, ax = plt.subplots(nrows = n_lipid_types, figsize = (figsize_x, 8*n_lipid_types), squeeze=False) 
	# the squeeze option is for when there is only one lipid type.  by default, subplots creates an axis object, not an array of axes, and so if only one lipid the indexing by lipid (lines 216 onwards) doesn't work. squeeze=False changes the default behaviour so axes can be indexed in all cases.
	# choose n ticks for x axis - this is ugly - could doubtless be improved but OK for now
	if nframes < 20:
		step = 4
	elif nframes < 50:
		step = 10
	elif nframes < 100:
		step = 20
	elif nframes < 400:
		step = 50
	elif nframes < 1000:
		step = 100
	elif nframes < 5000:
		step = 500
	elif nframes < 10000:
		step = 1000
	x_ticks_list = []
	for m in range(nmonomers_todisplay):
		x_ticks_list += range(m*nframes, (m+1)*nframes, step)
	for i, lipid in enumerate(lipid_list):
		ax[i,0].imshow(prot_dists_image[lipid], interpolation='nearest', cmap=cmx.RdPu, aspect='auto')
		ax[i,0].invert_yaxis()
		lipid_indices_to_show = lipid_indices_dict[lipid][lipids_to_show[lipid]]
		ax[i,0].set_yticks( numpy.arange(0, len(lipid_indices_to_show)*nsites, nsites) + 0.5 )
		ax[i,0].set_yticklabels(lipid_indices_to_show, rotation='horizontal')
		ax[i,0].set_ylabel('{} index'.format(lipid))
		ax[i,0].hlines(numpy.arange(0, len(lipid_indices_to_show)*nsites-1, nsites) + nsites - 0.45, -0.5,nframes*nmonomers_todisplay-0.5, colors='lightgray') # different from monomersnotsep
		plt.tick_params(axis='y', left=False, right=False)
		ax[i,0].set_xticks(numpy.array(x_ticks_list) - 0.5)
		ax[i,0].set_xticklabels(list(numpy.arange(0, nframes, step)*dt/1000)*nmonomers_todisplay)
		if i == (len(lipid_list) - 1):
			ax[i,0].set_xlabel('Time ($\mu$s)')
		ax[i,0].vlines(numpy.arange(0,nframes*nmonomers_todisplay, nframes) - 0.5, -0.5, len(lipid_indices_to_show)*nsites-0.5)
	plt.savefig(output_filename)
	#fig.close()

#def distances_image(prot_lipid_dists, nframes, lipid, nmonomers, lipid_indices, protein_index_list, sites_list, sites_index_colours_dict, dt=50, show_all=False, monomers_sep=False):
def distances_image(prot_lipid_dists_dict, nframes, lipid_list, nmonomers, lipid_indices_dict, protein_index_list, sites_list, sites_index_colours_dict, dt=50, show_all=False, monomers_sep=False):
	prot_lipid_dists_image_dict = {}
	lipid_string = lipid_list[0]
	for lipid in lipid_list[1:]:
		lipid_string += '_'+lipid
	for protein in protein_index_list:
		if monomers_sep:
			for m in range(nmonomers):
				for lipid in lipid_list:
					prot_lipid_dists_image_dict[lipid] = {}
					for site in sites_list:
						prot_lipid_dists_image_dict[lipid][site] = prot_lipid_dists_dict[lipid][protein][m][site]
				if show_all:
					output_filename = 'distance_analysis_{}_prot{}monomer{}_all_allres.png'.format(lipid_string, protein,m) # different from monomersnotsep
				else:
					output_filename = 'distance_analysis_{}_prot{}monomer{}_longints_allres.png'.format(lipid_string, protein,m) # different from monomersnotsep
				distances_image_single_graph(prot_lipid_dists_image_dict, nframes, lipid_list, nmonomers, lipid_indices_dict, sites_list, sites_index_colours_dict, output_filename, dt, show_all, monomers_sep)
		else:
			if show_all:
				output_filename = 'distance_analysis_{}_prot{}_all_allres.png'.format(lipid_string, protein)
			else:
				output_filename = 'distance_analysis_{}_prot{}_longints_allres.png'.format(lipid_string, protein)
			for lipid in lipid_list:
				prot_lipid_dists_image_dict[lipid] = {}
				for site in sites_list:
					prot_lipid_dists_image_dict[lipid][site] = numpy.concatenate(([prot_lipid_dists_dict[lipid][protein][m][site] for m in range(nmonomers)]), axis = 1)
			distances_image_single_graph(prot_lipid_dists_image_dict, nframes, lipid_list, nmonomers, lipid_indices_dict, sites_list, sites_index_colours_dict, output_filename, dt, show_all, monomers_sep)

def main(gro_file, trr_file, protein_name, lipid_list, prot_index_string, nframes, measure, stride=1, lipid_part='headgroup', cutoff=65, dt=50, show_all=False, monomers_sep=False, Nterm_start=0, compare='cutoff', structure='structure.gro', first_resnum=41):
	protein_seq = list(protein_info[protein_name]['protein_seq'])
	nmonomers = protein_info[protein_name]['protein_nmonomers']
	prot_index_list = make_prot_index_list(prot_index_string)
	#sites_list = protein_sites[protein_name].keys()
	#print sites_list
	#sites_index_colours_dict = sites_index_colours[protein_name]
	prot_lipid_dists_dict = {}
	lipid_indices_dict = {}
	for lipid in lipid_list:
		prot_lipid_dists_dict[lipid], lipid_indices_dict[lipid] = fetch_interactions(gro_file, trr_file, protein_seq, protein_name, nmonomers, lipid, prot_index_list, nframes, stride, lipid_part, cutoff)
		cov_matrix, close_cov_matrix, near_residues_all_frames, residue_interaction_strength_count = interaction_covariances(prot_lipid_dists_dict[lipid], protein_seq, nmonomers, len(lipid_indices_dict[lipid]), nframes, measure, compare)
		covariance_or_correlation_matrix_image(cov_matrix, measure, compare)
		covariance_or_correlation_matrix_image_closeonly(close_cov_matrix, near_residues_all_frames, Nterm_start, measure, compare)
		dendrogram_image(close_cov_matrix, near_residues_all_frames, Nterm_start, compare)
		network_based_on_correlation_analysis(cov_matrix, residue_interaction_strength_count, first_resnum=41)
		#clusters_on_structure(structure, cov_matrix, near_residues_all_frames, Nterm_start, compare)
		#distances_image(prot_lipid_dists_dict, nframes, lipid_list, nmonomers, lipid_indices_dict, prot_index_list, sites_list, sites_index_colours_dict, dt, show_all, monomers_sep)
		#dendrogram_image(cov_matrix, near_residues_all_frames, Nterm_start, compare)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Calculate lipid-to-interaction-site-distance covariance', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = parser.add_subparsers()

	parser_ints = subparsers.add_parser('cov', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_ints.add_argument('-f',help='input gro file', required=True)
	parser_ints.add_argument('-x', help='List of xtc files', nargs='+', required=True)
	parser_ints.add_argument('--protlist',default='0', help='''protlist should be input in the format:  0:4,6:8 (no spaces, index ranges are inclusive) 
	Protein numbering starts at 0)''')
	parser_ints.add_argument('-p', help='Protein name', choices=protein_info.keys(), required=True)

	parser_ints.add_argument('-lipid', nargs='+', help='Name of lipid', choices=lipid_particles['headgroup'], required=True)
	parser_ints.add_argument('-nframes', type=int, default = 400, help='Number of frames in trajectory')
	parser_ints.add_argument('-stride', type=int, default = 1, help='Frame number intervals that xtc will be read at')
	parser_ints.add_argument('-lp', default='headgroup', help='Part of lipid to consider in interactions', choices=lipid_particles.keys())
	parser_ints.add_argument('-d', type=float, default=6.5, help='Interactions cutoff distance (Angstroms)')
	parser_ints.add_argument('-dt', type=int, default=50, help='Number of nanoseconds per frame')
	parser_ints.add_argument('--showall', action='store_true', help='specificy to show all lipids, rather than just those that have an interaction with the protein')
	parser_ints.add_argument('--monomers', action='store_true', help='produce a separate graph for each monomer')
	parser_ints.add_argument('--cov', action='store_true', help='Produce a graph of covariance, rather than correlation coefficient (the default)')
	parser_ints.add_argument('--nterm', type=int, default=0, help='Integer value to add on residue numbers, in order to correct for residue number always starting at 1 in martini')
	#parser_ints.set_defaults(func=get_ints)
	parser_ints.add_argument('-compare', default='cutoff', help='correlation calculated either using binary on/off based on a cutoff dist, or based on distances between Chol and each residue', choices=['cutoff', 'dists'])
	parser_ints.add_argument('-s', help='structure to use to map clusters onto (gro file)')
	parser_ints.add_argument('-r', type=int, default = 41, help='Official residue number of the first residue in your structure')
	
	args = parser.parse_args()
	f = open('inputs.log', 'w')
	f.write('Args used: {}'.format(args)+'\n'+'command used: '+ ' '.join(sys.argv))
	f.close()

	if args.cov:
		measure = 'covariance'
	else:
		measure = 'correlation'
		
	main(args.f, args.x, args.p, args.lipid, args.protlist, args.nframes, measure, args.stride, args.lp, args.d, args.dt, args.showall, args.monomers, args.nterm, args.compare, args.s)

#	parser.add_argument('-i','--interaction_files', help='List of interaction files', nargs='+', type = lambda s: [item for item in s.split(' ')])
#	parser.add_argument('-s','--systems_list', help='List of system names given in interaction files list', nargs='+', type = lambda s: [item for item in s.split(' ')]


