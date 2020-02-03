'''Aug 2019 - added code that would allow for monomers to be treated separately'''

import os
import re
import ast
import numpy
import sys
import argparse

from interaction_freq_core import protein_info



files_dict = {
	'PM':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPI_PM_3x3_ARCHER/processed_trajectories/system_3x3_fromem1x1.noW.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPI_PM_3x3_ARCHER/processed_trajectories/md_prod3x3.until20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPI_PM_3x3_ARCHER/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'nprot':9, 'nframes':50000, 'dt':1},
	'P2subPE':{
		'GRO':'/sansom/n42/bioc1280/Kir2_2_3SPI_PIP2subPE_3x3/processed_trajectories/md_prod3x3.noW.firstframe.gro',
		'XTC':'/sansom/n42/bioc1280/Kir2_2_3SPI_PIP2subPE_3x3/processed_trajectories/md_prod3x3.until10us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir2_2_3SPI_PIP2subPE_3x3/processed_trajectories/md_prod3x3.from10us.to20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir2_2_3SPI_PIP2subPE_3x3/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'nprot':9, 'nframes':50000, 'dt':1},
	'PMnoGM3':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noGm3_3x3/processed_trajectories/md_prod3x3.noW.firstframe.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noGm3_3x3/processed_trajectories/md_prod3x3.until10us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPIsys3_noGm3_3x3/processed_trajectories//md_prod3x3.from10us.to20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPIsys3_noGm3_3x3/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'nprot':9, 'nframes':50000, 'dt':1},
	'PMnoCHOL':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noCHOL3_3x3/processed_trajectories/md_prod3x3.until20us.noW.firstframe.pbcmol.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noCHOL3_3x3/processed_trajectories/md_prod3x3.until20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPIsys3_noCHOL3_3x3/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'nprot':9, 'nframes':50000, 'dt':1},
	'PMnoP2GM':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPI_noP2noGM3_3x3_ARCHER/processed_trajectories/md_prod3x3.until20us.noW.firstframe.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPI_noP2noGM3_3x3_ARCHER/processed_trajectories/md_prod3x3.until20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPI_noP2noGM3_3x3_ARCHER/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'nprot':9, 'nframes':50000, 'dt':1},
	'PMnoP2GMCh':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noP2noGm3noCHOL2_3x3/processed_trajectories/md_prod3x3.until20us.noW.firstframe.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPIsys3_noP2noGm3noCHOL2_3x3/processed_trajectories/md_prod3x3.until20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPIsys3_noP2noGm3noCHOL2_3x3/processed_trajectories/md_prod3x3.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'nprot':9, 'nframes':50000, 'dt':1},
	'P2noPS':{
		'GRO':'/sansom/n42/bioc1280/Kir22_3SPIsys3_P2noPS_3x3/processed_trajectories/md_prod3x3_P2noPS.from0us.to20us.noW.firstframe.gro',
		'XTC':'/sansom/n42/bioc1280/Kir22_3SPIsys3_P2noPS_3x3/processed_trajectories/md_prod3x3_P2noPS.from0us.to20us.noW.dt1ns.pbcmol.center.xtc /sansom/n42/bioc1280/Kir22_3SPIsys3_P2noPS_3x3/processed_trajectories/md_prod3x3_P2noPS.from20us.to50us.noW.dt1ns.pbcmol.center.xtc',
		'nprot':9, 'nframes':50000, 'dt':1},
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
	

monomer_format_dict = { 
	'Kir2.2_3spi': {0:{'L':3, 'R':2}, 1:{'L':2, 'R':3}, 2:{'L':0, 'R':1}, 3:{'L':1, 'R':0}}
	}
# monomer format dict lists, for each monomer on the protein, which monomer is adjacent on the left and which is adjacent on the right
# this is so that modularity networks which have been produced giving only resids of the monomer can be annotated to indicate which interactions actually take place on adjacent monmers - and then resids for the whole protein structure can be used in calculations for residence times and occupancies, in order that only one site at a time is considered
# so the 'left' and 'right' tags correspond to the type of annotation you decide to add to modularity_network_annotated_ files
# the assumption is that each monomer only borders two other (as is the case for potassium channels)
# also assumed is that the protein monomers have incresing residue numbers per monomer - ie. residue numbers do not start at 1 again for each new protein monomer in the protein
	

def main(permonomer, nres_monomer, n_monomers, gro_filename, xtc_filename, nprot, lipid, nframes, dt, stride, lipid_part, monomer_format, d):
	filenames = os.listdir('.')
	for filename in filenames:
		if permonomer:
			matchobject = r'modularity_network_annotated_(\d+)\.txt'
		else:
			matchobject = r'modularity_network_(\d+)\.txt'
		if re.match(matchobject, filename):
			print filename
			annotation = '_'+re.match(matchobject, filename).group(1)
			f = open(filename)
			line_of_interest = f.readlines()[1]
			reslist = numpy.array(ast.literal_eval(line_of_interest))
			reslist_whole_lists = []
			if permonomer:
				for monomer in range(n_monomers):
					whole_list = []
					for res in reslist:
						resid = int(res[:-1])
						tag = res[-1]
						if tag == 'c':
							whole_list.append(resid + nres_monomer*monomer)
						elif tag == 'R':
							whole_list.append(resid + nres_monomer*monomer_format[monomer]['R'] )
						elif tag == 'L':
							whole_list.append(resid + nres_monomer*monomer_format[monomer]['L'] )
					#whole_list += list(reslist + nres_monomer*monomer)
					reslist_whole_lists.append(whole_list)
			else:
				whole_list = []
				for monomer in range(n_monomers):
					whole_list += list(reslist + nres_monomer*monomer)
				reslist_whole_lists.append(whole_list)
			for m,reslist in enumerate(reslist_whole_lists):
				if permonomer:
					annotation += '_m{}'.format(m)
				print reslist
				reslist_input = ''
				for res in whole_list:
					reslist_input += str(res)+','
				reslist_input = reslist_input[:-1]
				print reslist_input, annotation
				os.system('python /sansom/s105/bioc1280/Simulations/Tools/residence_time_lipid/lipid_res_time_v2.0.py get_ints -f {} -x {} -c reslist -reslist {} -nres {} -r {} -lipid {} -nframes {} -dt {} -stride {} --lipid_part headgroup -d {} -o {}'.format(gro_filename,xtc_filename, reslist_input, nres_monomer*n_monomers, nprot, lipid, nframes, dt, stride, d, annotation))
				# 6.5 used as cutoff distance since this is used in the distance covariance analysis also.

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Run restimes analysis on each modularity network', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('-sys', help='system name - used to fetch gro and xtc filenames, and nprot, dt and nframes', choices=files_dict.keys(), default='none')
	parser.add_argument('-p', help='protein name - used to fetch nresm and nmono', choices=protein_info.keys(), default='none')
	parser.add_argument('-f',help='input gro file - can be pre-specified in files_dict', default = 'none')
	parser.add_argument('-x', help='List of xtc files - can be pre-specified in files_dict', nargs='+', default = 'none')
	parser.add_argument('-np',default=1, help='Number of proteins in system - can be pre-specified in files_dict')
	parser.add_argument('-lipid', help='Name of lipid', required=True)
	parser.add_argument('-nframes', type=int, default = 50000, help='Number of frames in trajectory - can be pre-specified in files_dict')
	parser.add_argument('-stride', type=int, default = 1, help='Frame number intervals that xtc will be read at')
	parser.add_argument('-lp', default='headgroup', help='Part of lipid to consider in interactions: headgroup or phosphate')
	parser.add_argument('-d', type=float, default=6.5, help='Interactions cutoff distance (Angstroms)')
	parser.add_argument('-dt', type=int, default=1, help='Number of nanoseconds per frame (once "stride" has been taken into account) - can be pre-specified in files_dict')
	parser.add_argument('-permonomer', action='store_true', help='process results separately for each monomer')
	parser.add_argument('-nresm', type=int, default = 0, help='number of residues in each prot monomer - can be pre-specified in protein_info')
	parser.add_argument('-nmono', type=int, default = 0, help='Number of monomers per protein - can be pre-specified in protein_info')
	
	args = parser.parse_args()
	f = open('inputs_modularity_run_bs.log', 'w')
	f.write('Args used: {}'.format(args)+'\n'+'command used: python'+' '.join(sys.argv))
	f.close()

	# set nres_monomer and n_monomers - if numbers are given in command line this will override info in protein_info
	system = args.sys
	protname = args.p
	if args.permonomer and protname not in monomer_format_dict:
		print 'Add {} to monomer_format_dict'
	elif args.permonomer:
		monomer_format = monomer_format_dict[protname]
	else: 
		monomer_format = None
	if protname in protein_info and args.nresm == 0:
		nres_monomer = len(protein_info[protname]['protein_seq'])
	else:
		nres_monomer = args.nresm
	if protname in protein_info and args.nresm == 0:
		n_monomers = protein_info[protname]['protein_nmonomers']
	else: 
		n_monomers = args.nmonomer
	if args.f == 'none' and system not in files_dict:
		print 'Add gro filename to files_dict or specify using the -f flag'
	elif args.f == 'none':
		gro_filename = files_dict[system]['GRO']
	else:
		gro_filename = args.f
	if args.x == 'none' and system not in files_dict:
		print 'Add xtc filename(s) to files_dict or specify using the -x flag'
	elif args.x == 'none':
		xtc_filename = files_dict[system]['XTC']
		nprot = files_dict[system]['nprot']
		nframes = files_dict[system]['nframes']/args.stride
		dt = files_dict[system]['dt']/args.stride
	else:
		xtc_filename = ' '.join(args.x)
		nprot = args.np
		nframes = args.nframes
		dt = args.dt

	main(args.permonomer, nres_monomer, n_monomers, gro_filename, xtc_filename, nprot, args.lipid, nframes, dt, args.stride, args.lp, monomer_format, args.d)
