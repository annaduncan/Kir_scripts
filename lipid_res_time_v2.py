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
import os
import re
import numpy
import time
import MDAnalysis
from MDAnalysis.core.parallel.distances import distance_array
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cmx
import matplotlib.colors as colors

from scipy.optimize import curve_fit

import time

import ast
import argparse

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
lipid_colours = {'PIP2':'orange',\
                'PI3':'red',\
                'GM3':'magenta',\
                'CHOL':'cyan',\
                'PPCS':'darkgreen',\
                'POPE':'lightgrey',\
                'POPC':'black',\
                'POPS':'blue', 'MLCA':'green','CDL2':'yellow','CDL4':'yellow',
                'protein':'pink'}

def centroid_array_production_protein(protein_sel,num_protein_copies):
    dictionary_centroid_arrays = {}
    full_protein_coord_array = protein_sel.coordinates()
    list_individual_protein_coord_arrays = numpy.split(full_protein_coord_array,num_protein_copies)
    list_per_protein_centroids = [numpy.average(protein_coord_array,axis=0) for protein_coord_array in list_individual_protein_coord_arrays]
    dictionary_centroid_arrays['protein'] = numpy.array(list_per_protein_centroids)
    return dictionary_centroid_arrays


### main function to count interactions and produce interactions file ###
def count_frequencies(gro_file, trr_file, protein_res_total, protein_residue_list, nrepeats, lipid, nframes, dt, stride=1, lipid_part='headgroup', protein_centre='centroid', cutoff=6.5, annotation=''):

	universe = MDAnalysis.Universe(gro_file, trr_file)

	lipid_selection = 'resname {} and ('.format(lipid)
	for bead in lipid_particles[lipid_part][lipid]:
		lipid_selection += 'name {} or '.format(bead)
	lipid_selection = lipid_selection[:-4] + ')'
	
	lipid_rep_selection = 'resname {} and name {}'.format(lipid, lipid_particles[lipid_part][lipid][0]) # ie. just choose one bead

	lipids = universe.selectAtoms(lipid_selection)
	n_lipid_beads = len(lipid_particles[lipid_part][lipid])
	n_lipids = lipids.numberOfAtoms() / n_lipid_beads
	lipid_reps = universe.selectAtoms(lipid_rep_selection)
	#print 'lipid_reps.coordinates()', lipid_reps.coordinates()

    #initialise protein-lipid interactions list
	protein_lipid_interactions = numpy.empty( shape=(nrepeats,nframes), dtype=list )
	protein_lipid_restimes = dict(zip(range(nrepeats),[{} for i in range(nrepeats)]))
	print 'protein_lipid_restimes:',protein_lipid_restimes
	
	prot_lipid_occupancy = numpy.empty( shape=(nrepeats,nframes), dtype=list )
	
	startTime = time.time()    
	print 'Here we go...'
	frame = 0
	for ts in universe.trajectory[::stride]:
		if frame >= nframes:
			print 'Have reached maximum number of frames specified.  Stopping...'
			continue
		for i in range(nrepeats):
			single_prot = universe.segments[0][range(i*protein_res_total,(i+1)*protein_res_total)]
			# find protein centroid - or pick out residue to represent protein position
			if protein_centre == 'centroid':
				single_prot_cent = numpy.array([single_prot.centroid()])
			elif protein_centre == 'reslist':
				repeat_res_list = i*protein_res_total + numpy.array(protein_residue_list)
				single_prot_cent = universe.segments[0][repeat_res_list].coordinates()
			elif type(protein_centre) == int:
				#pick out BB of specified residue
				single_prot_cent = universe.segments[0][i*protein_res_total + protein_centre-1][0] # -1 is because res numbers are zero-indexed
				#print 'single_prot_cent', single_prot_cent
			else:
				print 'Error: protein_centre should either be an integer residue number, or "centroid"'
				return None
			if protein_centre == 'reslist':
				close = numpy.any(distance_array(single_prot_cent, lipid_reps.coordinates(), ts.dimensions) < cutoff, axis=0)
				nonzero_close = numpy.nonzero(close)[0]
				occupancy = numpy.sum(distance_array(single_prot_cent, lipid_reps.coordinates(), ts.dimensions) < cutoff, axis=0) > len(protein_residue_list)/2 ##does this work for proteins that are tetramers?
				prot_lipid_occupancy[i,frame] = np.any(occupancy)
			else:
				close = distance_array(single_prot_cent, lipid_reps.coordinates(), ts.dimensions) < cutoff
				nonzero_close = numpy.nonzero(close)[1] # numpy.nonzero gives a tuple of arrays - the second gives the indices of lipid residues that are 'close' to the protein
			for index_lipid in nonzero_close: 
				if protein_lipid_interactions[i,frame]:
					protein_lipid_interactions[i,frame].append(index_lipid)
				else:
					protein_lipid_interactions[i,frame] = [index_lipid]
				if index_lipid not in protein_lipid_restimes[i]:
					protein_lipid_restimes[i][index_lipid] = [dt]
				else:
					if frame>0 and protein_lipid_interactions[i,frame-1] and index_lipid in protein_lipid_interactions[i,frame-1]:
						protein_lipid_restimes[i][index_lipid][-1] += dt
					else:
						protein_lipid_restimes[i][index_lipid].append(dt)
#				if index_lipid not in protein_lipid_all_restimes_perlipid:
#					protein_lipid_all_restimes_perlipid[index_lipid] = [dt]
#				else:
#					if frame>0 and protein_lipid_interactions[i,frame-1] and index_lipid in protein_lipid_interactions[i,frame-1]:
#						protein_lipid_all_restimes_perlipid[index_lipid][-1] += dt
			#print i,protein_lipid_interactions[i,frame]
			#print i,protein_lipid_restimes[i]
			#update = '\rProtein {}/{} took {:3f} s'.format(i, nrepeats, time.time()-startTime)
			#print update,
			#sys.stdout.flush()
			#sys.stdout.write(update)
			#startTime = time.time()
		print 'Frame {} (fromtraj)or{} (hardcoded) took {:3f} s\r'.format(ts.frame, frame, time.time()-startTime),
		frame += 1
		startTime = time.time()
		#print protein_lipid_interactions[0:10,0:10], protein_lipid_restimes[0]
	print protein_lipid_interactions
	print protein_lipid_restimes
	protein_lipid_all_restimes_perprot = {}
	for prot,prot_list in protein_lipid_restimes.iteritems():
		all_lips = []
		for lipid_times in prot_list.values():
			all_lips += lipid_times
		protein_lipid_all_restimes_perprot[prot] = all_lips
	print protein_lipid_all_restimes_perprot 
	protein_lipid_all_restimes = []
	for protein_times in protein_lipid_all_restimes_perprot.values():
		protein_lipid_all_restimes += protein_times
	f = open('restime_data_{}{}'.format(lipid,annotation),'w')
	f.write(str(protein_lipid_interactions)+'\n')
	f.write(str(protein_lipid_restimes)+'\n')
	f.write(str(protein_lipid_all_restimes_perprot)+'\n')
	f.write(str(protein_lipid_all_restimes)+'\n')
	## save protein_lipid_occupancy
	f.close()
	return protein_lipid_interactions, protein_lipid_restimes, protein_lipid_all_restimes_perprot, protein_lipid_all_restimes


# functions for reading input - residue list and interactions file
def make_res_list(res_string, protein_res_total):
	# residue list should be input in the format 1:4,6:9,300,310:333
	# the following decodes that format into a straight list
	if res_string == 'all':
		res_list = numpy.arange(protein_res_total) + 1
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

#######################################
# Functions to make graphs of sigma and fit res time params:
#######################################

def read_restimes(lipid, annotation=''):
    f = open('restime_data_'+lipid+annotation)
    lines = f.readlines()
    f.close()
    return lines[-1].rstrip()
def read_all_restimes(lipids_of_interest, lipids_of_interest_annotations_dict={}):
	restimes = {}
	for lipid in lipids_of_interest:
		if lipid in lipids_of_interest_annotations_dict:
			for annotation in lipids_of_interest_annotations_dict[lipid]:
				restimes[lipid+annotation] = ast.literal_eval(read_restimes(lipid, annotation))
		else:
			restimes[lipid] = ast.literal_eval(read_restimes(lipid))
	return restimes

def find_nlips(top_file, lipids_of_interest):
    nlip = {}
    f = open(top_file)
    for line in f:
        for lipid in lipids_of_interest:
            if line[:len(lipid)] == lipid:
                nlip[lipid] = int(line.rstrip().split()[1])
    f.close()
    return nlip
    
def res_correlation_times(restimes, nframes, dt, lipids_of_interest_w_annotations, lipids_only_list, nlip):
    T = int(nframes*dt)
    delta_t_range = [0]+ range(5,100,5) + range(100,1000,50) + range(1000,T+1,100)
    sigma = {}
    for i,lipid in enumerate(lipids_of_interest_w_annotations):
        sigma[lipid] = {}
        restimes_j = numpy.array(restimes[lipid])
        sigma_0 = 1
        for delta_t in delta_t_range:
            if delta_t == 0:
                sigma[lipid][delta_t] = 1
                sigma0 = float(sum([restime - delta_t for restime in restimes_j[restimes_j>=delta_t]])) / ((T - delta_t) * nlip[lipids_only_list[i]])

            else:
                try:
                    sigma[lipid][delta_t] = float(sum([restime - delta_t for restime in restimes_j[restimes_j>=delta_t]])) / ((T - delta_t) * nlip[lipids_only_list[0]] * sigma0)
                except ZeroDivisionError:
                    sigma[lipid][delta_t] = 0
    return sigma   

def find_residence_time_tripleexp(time_array, sigma_array):
    if time_array.shape != sigma_array.shape:
        print 'Time and sigma arrays should have the same length'
        return None
    def model_function(t, A, B, C, d,f,g):
        return A*numpy.exp(-d*t) + B*numpy.exp(-f*t) + C*numpy.exp(-g*t)
    params, params_covariance = curve_fit(model_function, time_array, sigma_array, [1,0.1,1,1,1,1])
    sample_t = numpy.linspace(time_array[0], time_array[-1], 50)
    sample_sigma = model_function(sample_t, *params)
    
    param_std = numpy.sqrt(numpy.diagonal(params_covariance))
    first_decay_const = 1./(params[2]*1000.)
    second_decay_const = 1./(params[3]*1000.)
    third_decay_const  = 1./(params[4]*1000.)
    first_const_error = 1./(param_std[2]*1000000.)
    second_const_error = 1./(param_std[3]*1000000.) #1000 converts from ns to us
    third_const_error = 1./(param_std[4]*1000000.)
    
    return sample_t, sample_sigma, first_decay_const, first_const_error, second_decay_const, second_const_error, third_decay_const, third_const_error

def find_residence_time_doubleexp(time_array, sigma_array):
    if time_array.shape != sigma_array.shape:
        print 'Time and sigma arrays should have the same length'
        return None
    def model_function(t, A, B, c, d):
        return A*numpy.exp(-c*t) + B*numpy.exp(-d*t)
    params, params_covariance = curve_fit(model_function, time_array, sigma_array, [1,0.1,1,1])
    sample_t = numpy.linspace(time_array[0], time_array[-1], 50)
    sample_sigma = model_function(sample_t, *params)
    
    param_std = numpy.sqrt(numpy.diagonal(params_covariance))
    first_decay_const = 1./(params[2]*1000.)
    second_decay_const = 1./(params[3]*1000.) #1000 converts from ns to us
    first_const_error = 1./(param_std[2]*1000000.)
    second_const_error = 1./(param_std[3]*1000000.) #1000 000 converts from ns to us for std since this is a sqrt
    
    return sample_t, sample_sigma, first_decay_const, first_const_error, second_decay_const, second_const_error    

def find_residence_time_singleexp(time_array, sigma_array):
    if time_array.shape != sigma_array.shape:
        print 'Time and sigma arrays should have the same length'
        return None
    def model_function(t, A, b):
        return A*numpy.exp(-b*t)
    params, params_covariance = curve_fit(model_function, time_array, sigma_array, [1,1])
    sample_t = numpy.linspace(time_array[0], time_array[-1], 50)
    sample_sigma = model_function(sample_t, *params)
    
    param_std = numpy.sqrt(numpy.diagonal(params_covariance))
    first_decay_const = 1./(params[1]*1000.)
    first_const_error = 1./(param_std[1]*1000000.)
   
    return sample_t, sample_sigma, first_decay_const, first_const_error   

def get_residence_times(sigma, lipids_of_interest, nframes, dt, fittingfunc='doubleexp'):
    lipid_list = lipids_of_interest
    delta_t_range = numpy.array([0]+ range(5,100,5) + range(100,1000,50) + range(1000,int(nframes*dt+1),100),dtype=float)
    n = len(lipid_list)
    residence_times = {}
#    fig,axes = plt.subplots(1,n,figsize=(n*2.5,2),squeeze=False)
    for i,lipid in enumerate(lipid_list):
        residence_times[lipid] = {'sample_t':'', 'sample_sigma':'', 'first_decay_const':0., 'first_const_error':0., 
                                            'second_decay_const':0., 'second_const_error':0., 'third_decay_const':0., 'third_const_error':0.}
        sigma_array = numpy.array([sigma[lipid][delta_t] for delta_t in delta_t_range])
        if fittingfunc == 'singleexp':
            residence_times[lipid]['sample_t'], residence_times[lipid]['sample_sigma'], \
            residence_times[lipid]['first_decay_const'], residence_times[lipid]['first_const_error'] \
            = find_residence_time_singleexp(delta_t_range, sigma_array)
        elif fittingfunc == 'doubleexp':
            residence_times[lipid]['sample_t'], residence_times[lipid]['sample_sigma'], \
            residence_times[lipid]['first_decay_const'], residence_times[lipid]['first_const_error'], \
            residence_times[lipid]['second_decay_const'], residence_times[lipid]['second_const_error'] \
            = find_residence_time_doubleexp(delta_t_range, sigma_array)
        elif fittingfunc == 'tripleexp':
            residence_times[lipid]['sample_t'], residence_times[lipid]['sample_sigma'], \
            residence_times[lipid]['first_decay_const'], residence_times[lipid]['first_const_error'], \
            residence_times[lipid]['second_decay_const'], residence_times[lipid]['second_const_error'], \
            residence_times[lipid]['third_decay_const'], residence_times[lipid]['third_const_error'] \
            = find_residence_time_tripleexp(delta_t_range, sigma_array)
    return residence_times

def plot_restimes_and_realdata(sigma, residence_times, figname, lipids_of_interest_w_annotations, lipids_only_list, nframes, dt, fittingfunc='doubleexp'):
    lipid_list = lipids_of_interest_w_annotations
    delta_t_range = [0]+ range(5,100,5) + range(100,1000,50) + range(1000,int(nframes*dt+1),100)
    n = len(lipid_list)
    fig,axes = plt.subplots(n,1,figsize=(2.5,n*2),squeeze=False)
    for i,lipid in enumerate(lipid_list):
        axes[i,0].plot((numpy.array(delta_t_range,dtype=float)/1000.)[:-1], [sigma[lipid][delta_t] for delta_t in delta_t_range][:-1], 
                       color=lipid_colours[lipids_only_list[i]], linewidth=3)
        axes[i,0].plot((residence_times[lipid]['sample_t']/1000.), residence_times[lipid]['sample_sigma'],
                       color=lipid_colours[lipids_only_list[i]], linestyle = ':', linewidth=5)
        #plt.xlim((0,20050))
        axes[i,0].set_ylim((10**(-6),1))
        axes[i,0].set_yscale('log')
        #axes[i,0].set_title(lipid)
	if fittingfunc=='tripleexp':
		axes[i,0].text(5,0.001, r'{:5.2f} $\pm${:5.2f} $\mu s$'.format(residence_times[lipid]['third_decay_const'], residence_times[lipid]['third_const_error']))
	if fittingfunc in ['doubleexp','tripleexp']:
		axes[i,0].text(5,0.0001, r'{:5.2f} $\pm${:5.2f} $\mu s$'.format(residence_times[lipid]['second_decay_const'], residence_times[lipid]['second_const_error']))
	if fittingfunc in ['singleexp','doubleexp','tripleexp']:
		axes[i,0].text(5,0.00001, r'{:5.2f} $\pm${:5.2f} $\mu s$'.format(residence_times[lipid]['first_decay_const'], residence_times[lipid]['first_const_error']))
        #plt.grid(True,which="both",ls="-",color="lightgrey", zorder=-1)
        #ax.set_axisbelow(True)    
    #plt.show()
    plt.savefig(figname,bbox_inches='tight')
    #plt.close() 
    
def output_graphs_and_restimes(figname, txtfilename, top_file, lipids_of_interest, nframes, dt, fittingfunc_list=['doubleexp'], annotation_dict={}):
	restimes = read_all_restimes(lipids_of_interest, annotation_dict)
	nlip = find_nlips(top_file, lipids_of_interest)
	print 'nlip', nlip
	lipids_of_interest_w_annotations = []
	lipids_only_list = []
	for lipid in lipids_of_interest:
		if lipid in annotation_dict:
			for annotation in annotation_dict[lipid]:
				lipids_of_interest_w_annotations.append(lipid+annotation)
				lipids_only_list.append(lipid)
		else:
			lipids_of_interest_w_annotations.append(lipid)
			lipids_only_list.append(lipid)
	sigma = res_correlation_times(restimes, nframes, dt, lipids_of_interest_w_annotations, lipids_only_list, nlip)
	# order lipids_of_interest_w_annotations -to do
	# curve fitting
	for fittingfunc in fittingfunc_list:
		residence_times = get_residence_times(sigma, lipids_of_interest_w_annotations, nframes, dt, fittingfunc)
		figname_ff = figname[:-4] + fittingfunc + figname[-4:]
		txtfilename_ff = txtfilename[:-4] + fittingfunc + txtfilename[-4:]
		#plot_restimes_and_realdata(sigma, residence_times, figname_ff, lipids_of_interest_w_annotations, lipids_only_list, nframes, dt, fittingfunc)
		f = open(txtfilename_ff, 'w')
		for lipid in lipids_of_interest_w_annotations:
			if fittingfunc == 'singleexp':
				f.write('{}: {:5.2f} +/- {:5.2f} mu s\n'.format(lipid, residence_times[lipid]['first_decay_const'], residence_times[lipid]['first_const_error']))
			elif fittingfunc == 'doubleexp':
				f.write('{}: {:5.2f} +/- {:5.2f} mu s, {:5.2f} +/- {:5.2f} mu s\n'.format(lipid, residence_times[lipid]['first_decay_const'], residence_times[lipid]['first_const_error'], \
				residence_times[lipid]['second_decay_const'], residence_times[lipid]['second_const_error']))
			elif fittingfunc == 'tripleexp':
				f.write('{}: {:5.2f} +/- {:5.2f} mu s, {:5.2f} +/- {:5.2f} mu s, {:5.2f} +/- {:5.2f} mu s\n'.format(lipid, residence_times[lipid]['first_decay_const'], residence_times[lipid]['first_const_error'], \
				residence_times[lipid]['second_decay_const'], residence_times[lipid]['second_const_error'], \
				residence_times[lipid]['third_decay_const'], residence_times[lipid]['third_const_error']))
		f.close()



if __name__ == "__main__":
	# sub-command funcs for the parser - ie. workflows available when running script:
	def get_ints(args):
		f = open('get_ints.log', 'w')
		f.write('Args used: {}'.format(args)+'\n')
		f.close()
		res_list = make_res_list(args.reslist, args.nres)
		#protein_seq = list(protein_info[args.p]['protein_seq'])
		print args.c
		if args.c not in ['centroid', 'reslist']:
			args.c = int(args.c) # this should return an error if the argument is not integer-like
			if args.c < args.nres:
				#check that the residue number is valid
				print 'Error: residue number selected is larger than the number of residues in the protein sequence'
				return None
		for lipid in args.lipid:
			proteinres_lipid_interactions = count_frequencies(args.f, args.xtc, args.nres, res_list, args.nrepeats, lipid, args.nframes, args.dt, args.stride, args.lipid_part, args.c, args.cutoff, args.o)
	def plot_figs(args):
#		if args.figname == 'restimes.svg':
#			args.figname = 'restimes_'+args.fittingfunc+'.svg'
#		if args.txtfilename == 'final_restimes.txt':
#			args.txtfilename = 'final_restimes_'+args.fittingfunc+'.txt'
		if args.annotation != []:
			annotation_dict = {}
			annotation_dict[list(args.lipid)[0]] = list(args.annotation)
			print annotation_dict
		output_graphs_and_restimes(args.figname, args.txtfilename, args.top_file, list(args.lipid), args.nframes, args.dt, list(args.fittingfunc), annotation_dict)

		
	parser = argparse.ArgumentParser(description='Analyse prot-lipid interactions at the residue level', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = parser.add_subparsers()

	parser_ints = subparsers.add_parser('get_ints', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_ints.add_argument('-f',help='input gro file', required=True)
	parser_ints.add_argument('-x', '--xtc', help='List of xtc files', nargs='+', required=True)
	parser_ints.add_argument('-reslist',default='all', help='''reslist should be input in the format:  1:4,6:9,300,310:333  \n(no spaces, residue ranges are inclusive,  \n
	Residue numbering starts at 0 - eg. residue 1 in your gro file is residue 0 when inputting here.)\n
	OR ALL residues can be specified using: all''')
	parser_ints.add_argument('-nres', type= int, help='Number of protein residues', required=True)
	parser_ints.add_argument('-r', '--nrepeats', type=int, help='Number of protein MONOMERS in system', required=True)
	parser_ints.add_argument('-lipid', nargs='+', help='Name of lipid', choices=lipid_particles['headgroup'], required=True)
	parser_ints.add_argument('-nframes', type=int, default = 400, help='Number of frames in trajectory')
	parser_ints.add_argument('-dt', type=float, default = 5, help='Timestep of each frame')
	parser_ints.add_argument('-stride', type=int, default = 1, help='Frame number intervals that xtc will be read at')
	parser_ints.add_argument('-lp', '--lipid_part', default='headgroup', help='Part of lipid to consider in interactions', choices=lipid_particles.keys())
	parser_ints.add_argument('-c', default='centroid', help='Part of protein to measure lipid distances from, for first approximation - can be either "centroid" for protein centroid, or "reslist" if you are specifying a list of residues with -reslist, or an integer residue number (residue numbering starts from 1)')
	#parser_ints.add_argument('-cd', type = float, default=80, help='Distance away from protein to select lipids that will undergo a closer inspection (a lower number speeds the script up, but too low may mean that some lipids aren\'t considered')
	parser_ints.add_argument('-d', '--cutoff', type=float, default=45, help='Interactions cutoff distance (nm)')
	parser_ints.add_argument('-o', default='', help='Extra text to add to output filenames')
	parser_ints.set_defaults(func=get_ints)

	parser_plot = subparsers.add_parser('plot', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_plot.add_argument('-lipid', nargs='+', help='Name of lipid', choices=lipid_particles['headgroup'], required=True)
	parser_plot.add_argument('-nframes', type=int, default = 400, help='Number of frames in trajectory')
	parser_plot.add_argument('-annotation', nargs='+', default='', help='list of annotations - if using this option, there should only be a single lipid in the lipid list')
	parser_plot.add_argument('-dt', type=float, default = 5, help='Timestep of each frame')
	parser_plot.add_argument('-ofig', dest='figname', help='Name of output figure file', default='restimes.svg')
	parser_plot.add_argument('-odata', dest='txtfilename', help='Name of output data file', default = 'final_restimes.txt')
	parser_plot.add_argument('-top', dest='top_file', help='Name of top file that includes number of lipids in simulation', required=True)
	parser_plot.add_argument('-ff', nargs='+', dest='fittingfunc', default = 'doubleexp', help='Type of exponential function to fit data to.  Options are: singleexp, doubleexp or tripleexp. Possible to list more than one option')
	parser_plot.set_defaults(func=plot_figs)
	
	args = parser.parse_args()
	
	# record script usage:
	log_file_num = [0]
        for file in os.listdir('.'):
            m = re.match(r'commands_used.(\d+).log', file)
            if m:
                log_file_num.append(int(m.group(1)))
        logfile = 'commands_used{:03d}.log'.format(max(log_file_num)+1)
	f = open(logfile,'w')
	f.write('Command used:\npython '+(' ').join(sys.argv)+'\n')
	f.write('All arguments: \n'+str(args)+'\n')
	f.close()
	
	args.func(args)

#	parser.add_argument('-i','--interaction_files', help='List of interaction files', nargs='+', type = lambda s: [item for item in s.split(' ')])
#	parser.add_argument('-s','--systems_list', help='List of system names given in interaction files list', nargs='+', type = lambda s: [item for item in s.split(' ')]
