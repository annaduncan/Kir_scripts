import os
import re
import ast
import numpy
import sys

from interaction_freq_core import protein_info

def make_tcl_file(gro_filename, nres_monomer, n_monomers, system, lipid, annotations_list, fitting, permonomer):
	#colours = ['magenta', 'green', 'blue', 'red', 'orange', 'yellow', 'tan', 'pink', 'cyan' , 'cyan3', 'blue2', 'violet', 'iceblue', 'green3', 'cyan2', 'blue3', 'violet2', 'red2', 'red3', 'orange2', 'orange3', 'purple', 'lime']
	#colours = ['magenta', 'green', 'red', 'magenta', 'green', 'blue', 'red', 'orange', 'magenta', 'violet2', 'green', 'blue', 'red', 'orange','magenta', 'violet2','green', 'red', 'iceblue']
	colours = ['yellow', 'blue', 'red', 'tan', 'pink', 'cyan',  'iceblue', 'green', 'colour2', '6','8','11','12','13','14','16','17','18','19','20']
	#'orange', 'yellow', 'tan', 'pink', 'cyan' , 'cyan3', 'blue2', 'violet', 'iceblue', 'green3', 'cyan2', 'blue3', 'violet2', 'red2', 'red3', 'orange2', 'orange3', 'purple', 'lime']
	
	if annotations_list == []:
		filenames = os.listdir('.')
		annotations_list = []
		for filename in filenames:
			if re.match(r'modularity_network_(\d+)\.txt', filename):
				print filename
				annotation = re.match(r'modularity_network_(\d+)\.txt', filename).group(1)
				annotations_list.append(annotation)
	
		# make dictionary of residence times, then sort annotation files from longest to shortest residence times

		restimes_dict = {}
		f = open('final_restimes{}exp.txt'.format(fitting))
		lines = f.readlines()
		f.close()
		print 'lines: ',lines
		for line in lines:
			if fitting == 'single':
				if re.match(lipid + r'_(\d+): (\d*[\d\.\s]{5}) \+/\- \d*[\d\.\s]{5} mu s', line):
					m = re.match(lipid + r'_(\d+): (\d*[\d\.\s]{5}) \+/\- \d*[\d\.\s]{5} mu s', line)
			elif fitting == 'double':
				if re.match(lipid + r'_(\d+): \d*[\d\.\s]{5} \+/\- \d*[\d\.\s]{5} mu s, (\d*[\d\.\s]{5}) \+/\- \d*[\d\.\s]{5} mu s', line):
					m = re.match(lipid + r'_(\d+): \d*[\d\.\s]{5} \+/\- \d*[\d\.\s]{5} mu s, (\d*[\d\.\s]{5}) \+/\- \d*[\d\.\s]{5} mu s', line)
					print 'matched line: ', line
			elif fitting == 'triple':
				if re.match(lipid + r'_(\d+): \d*[\d\.\s]{5} \+/\- \d*[\d\.\s]{5} mu s, \d*[\d\.\s]{5} \+/\- \d*[\d\.\s]{5} mu s, (\d*[\d\.\s]{5}) \+/\- \d*[\d\.\s]{5} mu s', line):
					m = re.match(lipid + r'_(\d+): \d*[\d\.\s]{5} \+/\- \d*[\d\.\s]{5} mu s, \d*[\d\.\s]{5} \+/\- \d*[\d\.\s]{5} mu s, (\d*[\d\.\s]{5}) \+/\- \d*[\d\.\s]{5} mu s', line)
			else:
				print '3rd argument should be one of: single, double, triple'
			restimes_dict[m.group(1)] = float(m.group(2))
		annotations_by_restime = sorted(restimes_dict.items(), key = lambda kv:(kv[1],kv[0]),reverse=True)
		f = open('restimes_on_struct_restimes_and_colours_{}_{}_{}fitting.txt'.format(system, lipid, fitting), 'w')
		for i,a in enumerate(annotations_by_restime):
			f.write(str(a)+ ' '+ colours[i] + '\n')
		f.close()
	else:
		annotations_by_restime = zip(annotations_list, [0]*len(annotations_list))

	#vmd_color_list = [27,7,0,1,3,4,5,9,10,22,23,25,15,20,21,24,26,29,30,31,32,11,12]
	#vmd_color_list = [27,7,1,27,7,0,1,3,27,25,7,0,1,3,27,25,7,1,10]
	vmd_color_list = [4,0,1,5,9,10,15,7,2,6,8,11,12,13,14,16,17,18,19,20]
	#3,4,5,9,10,22,23,25,15,20,21,24,26,29,30,31,32,11,12]
# colours: magenta, green, blue, red, orange, yellow, tan, pink cyan , cyan3, blue2, violet, iceblue, green3, cyan2,...purple, lime, 
	f = open('view_binding_sites_{}_{}.vmd'.format(system,lipid),'w')
	f.write('display set view 600 900\n')
	f.write('mol addrep 0\n')
	f.write('mol new {{{}}} type {{gro}} first 0 last -1 step 1 waitfor 1\n'.format(gro_filename))
	f.write('display rendermode GLSL\n')
	f.write('axes location Off\n')
	f.write('display depthcue off\n')
	f.write('material change outline AOChalky 0.390000\n')
	f.write('material change outlinewidth AOChalky 0.110000\n')
	f.write('material change outline AOChalky 1.020000\n')
	f.write('material change outline AOChalky 2.070000\n')
	f.write('material change outlinewidth AOChalky 0.160000\n')
	f.write('mol material AOChalky\n')
	f.write('mol representation QuickSurf 1.100000 1.500000 0.500000 1.000000\n')
	f.write('rotate x by -90.000000\n')
	f.write('mol modselect 0 0 name BB\n')
	f.write('mol modstyle 0 0 QuickSurf 1.300000 1.500000 0.500000 1.000000\n')
	f.write('mol modmaterial 0 0 AOChalky\n')
	f.write('mol modcolor 0 0 ColorID 8\n')
	f.write('mol representation QuickSurf 1.400000 1.500000 0.500000 1.000000\n')
	for i,annotationtuple in enumerate(annotations_by_restime):
		annotation = annotationtuple[0]
		if permonomer:
			af = open('modularity_network_monomersep_{}.txt'.format(annotation))
		else:
			af = open('modularity_network_{}.txt'.format(annotation))
		line_of_interest = af.readlines()[1]
		reslist = numpy.array(ast.literal_eval(line_of_interest)) + 1 # +1 is because MDanalysis reads 1st residue as residue 0, whereas gro file has 1st residue as residue 1
		whole_list = []
		for monomer in range(n_monomers):
			whole_list += list(reslist + nres_monomer*monomer)
		print reslist
		reslist_input = ''
		for res in whole_list:
			reslist_input += str(res)+' '
		reslist_input = reslist_input[:-1]
		print reslist_input, annotation
		f.write('mol addrep 0\n')
		f.write('mol modselect {} 0 resid {}\n'.format(i+1,reslist_input))
		f.write('mol modcolor {} 0 ColorID {}\n'.format(i+1, vmd_color_list[i]))
	f.write('set viewpoints([molinfo top]) {{{1 0 0 -57.5827} {0 1 0 -57.3361} {0 0 1 -101.682} {0 0 0 1}} {{-0.708475 0.70564 0.0114819 0} {0.00638021 -0.0098648 0.999931 0} {0.705705 0.708499 0.0024869 0} {0 0 0 1}} {{0.0219772 0 0 0} {0 0.0219772 0 0} {0 0 0.0219772 0} {0 0 0 1}} {{1 0 0 0.04} {0 1 0 0.15} {0 0 1 0} {0 0 0 1}}}\n')
#	set viewpoints([molinfo top]) {{{1 0 0 -57.5827} {0 1 0 -57.3361} {0 0 1 -101.682} {0 0 0 1}} {{-0.708475 0.70564 0.0114819 0} {0.00638021 -0.0098648 0.999931 0} {0.705705 0.708499 0.0024869 0} {0 0 0 1}} {{0.0219772 0 0 0} {0 0.0219772 0 0} {0 0 0.0219772 0} {0 0 0 1}} {{1 0 0 0.04} {0 1 0 0.15} {0 0 1 0} {0 0 0 1}}}
	f.write('lappend viewplist [molinfo top]\n')
	f.write('set topmol [molinfo top]\n')
	f.write('# done with molecule 0\n')
	f.write('foreach v $viewplist {\n')
	f.write('  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)\n')
	f.write('}\n')
	f.write('foreach v $fixedlist {\n')
	f.write('  molinfo $v set fixed 1\n')
	f.write('}\n')
	f.write('unset viewplist\n')
	f.write('unset fixedlist\n')
	f.write('mol top $topmol\n')
	f.write('unset topmol\n')

	f.write('render Tachyon binding_sites_{}_{} "/sbcb/packages/opt/Linux_x86_64/vmd/1.9.2/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -o %s.tga'.format(system,lipid))
	f.close()

#nres_monomer=332
#n_monomers=4
#nres_monomer=1328
#n_monomers=1
#gro_filename= '/sansom/s105/bioc1280/Simulations/Kir2_2/sys3_3SPI/exchange_lipids1/test_em_1x1_10ns.gro'

print 'Use: system, lipid, protein, gro_filename, permonomer, fitting  <list of modularity network file indexes>'
print 'protein can be one of: {}'.format(protein_info.keys())

system = sys.argv[1]
lipid = sys.argv[2]
protein = sys.argv[3]
gro_filename = sys.argv[4]
permonomer = bool(int(sys.argv[5]))
fitting = sys.argv[6]
#option to supply annotations as list of numbers if already know the order you want these to load in
if len(sys.argv) > 7:
	annotations_list = sys.argv[7:]
else:
	annotations_list = []

nres_monomer = len(protein_info[protein]['protein_seq'])
n_monomers = protein_info[protein]['protein_nmonomers']

print 'args used: {}'.format(sys.argv)

make_tcl_file(gro_filename, nres_monomer, n_monomers, system, lipid, annotations_list, fitting, permonomer)
