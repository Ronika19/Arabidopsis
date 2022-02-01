import re

file1 = open('WGCNA_Modules.txt','r')
line1 = file1.readline()
line1 = file1.readline()
genes, modules = [],[];
while line1:
	line1 = line1.rstrip()
	split_line1 = line1.split()
	genes.append(split_line1[0])
	modules.append(int(split_line1[1]))
	line1 = file1.readline()
print(len(genes),len(modules))

mods = []
for module in modules:
	if module not in mods:
		mods.append(int(module))
mods.sort(); print(mods);
x = 0; gene=[];

file2 = open('Arabidopsis_Modules_Clusters.txt','w')
file3 = open('WGCNA_Arabidopsis_Modules_Clusters.txt','w')
for mod in mods:
	print(mod); file2.write(str(mod)+"\t")
	for x in range(len(modules)):
		if (modules[x] == mod):
			print(mod, modules[x], genes[x])
			gene.append(genes[x])
			file2.write(genes[x]+"\t")
			file3.write(genes[x]+"\t")
			#file2.write(str(mod)+"\t"+str(modules[x])+"\t"+genes[x]+"\n")
	file2.write("\n")
	file3.write("\n")

