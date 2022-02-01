import glob

# All DEGs
f1 = glob.glob('Arabidopsis_Cluster_DEG_*.txt')
#file1 = open('Arabidopsis_Cluster_DEG_Comp1_CTvsHL.txt','r')
file1 = open(f1[0], 'r')
line1 = file1.readline()
DEG_transcripts = []; indices = [];
while line1:
	line1 = line1.rstrip()
	split_line1 = line1.split('\t')
	for tids in split_line1:
		DEG_transcripts.append(tids)
	DEG_transcripts.append('\n')
	line1 = file1.readline()
print(len(DEG_transcripts))

for items in range(len(DEG_transcripts)):
	if DEG_transcripts[items] == '\n':
		indices.append(items); print(items)

# Modules using WGCNA
file2 = open('WGCNA_Modules.txt','r')
for i in range(2):
	line2 = file2.readline()
modules, gene_id = [],[];
while line2:
	line2 = line2.rstrip()
	split_line2 = line2.split('\t')
	gene_id.append(split_line2[0].replace('"',''))
	modules.append(split_line2[1])
	line2 = file2.readline()
#print(DEG_transcripts,"\n",gene_id)

file3 = open('Arabidopsis_DEG_Modules.txt','w')
#file4 = open('DC_DEG_Modules.txt','w')
#file5 = open('LC_DEG_Modules.txt','w')
k = 0; DEG1_Module = [];
for genes in DEG_transcripts:
	if (k < int(indices[0])):	# Cluster of Pathogenesis DEGs
		if (genes in gene_id):
			indexes = gene_id.index(genes); #print(indexes, genes, gene_id[indexes]);
			DEG1_Module.append(modules[indexes])
			file3.write(gene_id[indexes]+"\t"+modules[indexes]+"\n")
		#if (genes not in gene_id): # DEGs not in Modules
		#	print(indexes, genes);
	k += 1

#print(DEG1_Module,"\n\n", DEG2_Module,"\n\n", DEG3_Module)
file3.close(); file1.close(); file2.close();

file4 = open('Arabidopsis_DEG_Modules.txt','r')
line4 = file4.readline()
geneid, mod = [],[];
while line4:
	line4 = line4.rstrip()
	split_line4 = line4.split('\t')
	geneid.append(split_line4[0])
	mod.append(int(split_line4[1]))
	line4 = file4.readline()

modset = sorted(set(mod), reverse=False); #print(modset);

file5 = open('Arabidopsis_DEG_Genes_Modules.txt','w')
for i in modset:
	file5.write(str(i)+'\t')
	indices = [index for index, element in enumerate(mod) if element == i]; #print(i, indices);
	m = 0
	for j in indices:
		m += 1; #print(geneid[j], mod[j])
		if m < len(indices):
			file5.write(geneid[j]+',')
		elif m == len(indices):
			file5.write(geneid[j]+'\n')



