import re
import glob

infiles = ['WGCNA_Modules.txt', 'Comp1_CTvsHL.txt']
outfiles = ['Arabidopsis_Modules_Clusters.txt', 'WGCNA_Arabidopsis_Modules_Clusters.txt', 'Arabidopsis_Cluster_DEG_Comp1_CTvsHL.txt']

def Module_Cluster():
	file1 = open(infiles[0],'r')
	line1 = file1.readline()
	line1 = file1.readline()
	genes, modules = [],[];
	while line1:
		line1 = line1.rstrip()
		split_line1 = line1.split()
		genes.append(split_line1[0])
		modules.append(int(split_line1[1]))
		line1 = file1.readline()

	mods = []
	for module in modules:
		if module not in mods:
			mods.append(int(module))
	mods.sort(); 
	x = 0; gene=[];

	file2 = open(outfiles[0],'w')
	file3 = open(outfiles[1],'w')
	for mod in mods:
		file2.write(str(mod)+"\t")
		for x in range(len(modules)):
			if (modules[x] == mod):
				gene.append(genes[x])
				file2.write(genes[x]+"\t")
				file3.write(genes[x]+"\t")
		file2.write("\n")
		file3.write("\n")

def DEG_Cluster():
	file1 = open(infiles[1],'r')
	file4 = open(outfiles[2],'w')

	for i in range(2):
		line1 = file1.readline()
	while line1:
		line1 = line1.rstrip()
		split_line1 = line1.split('\t')
		genes1 = split_line1[0]
		fold_change1 = split_line1[2]
		if (split_line1[6] != 'NA'):
			if (float(split_line1[6]) < 0.05):
				file4.write(genes1+"\t")
		line1 = file1.readline()
	file4.write("\n")
	file4.close()





