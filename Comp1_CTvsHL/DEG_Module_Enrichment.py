import scipy.stats as stats

nodes_file = open('Arabidopsis_Nodes.txt','r')
line = nodes_file.readlines()
len_nodes = len(line[1:]); #print(len_nodes);	# Total no of Genes in WGCNA Network

file1_1 = open('WGCNA_Arabidopsis_Modules_Clusters.txt','r')
Total_Modules = len(file1_1.readlines()); #print(Total_Modules)
file1_1.close()
file1 = open('WGCNA_Arabidopsis_Modules_Clusters.txt','r')
line1 = file1.readline()
Mod_Size, Mod_Size_Percent = [],[]
while line1:
	line1 = line1.rstrip()
	split_line1 = line1.split('\t'); #print(len(split_line1), (len(split_line1)/len_nodes)*100);
	Mod_Size.append(len(split_line1))	# No of Genes or Gene Count in each WGCNA Modules = Mod_Size[1:]
	Mod_Size_Percent.append((len(split_line1)/len_nodes)*100)	# Percentage of Genes in each WGCNA Modules = Mod_Size_Percent[1:]
	line1 = file1.readline()

file2 = open('Arabidopsis_DEG_Modules.txt','r')
line2 = file2.readline()
Arabidopsis_dict = {}; len_Arabidopsis_degs = 0
while line2:
	line2 = line2.rstrip()
	split_line2 = line2.split('\t')
	Arabidopsis_deg_mods = int(split_line2[1]); #print(Arabidopsis_deg_mods);
	for i in range(Total_Modules):
		if (Arabidopsis_deg_mods == i):
			Arabidopsis_dict[split_line2[0]] = i
	len_Arabidopsis_degs += 1	# No of Arabidopsis DEGs
	line2 = file2.readline()
#print(Arabidopsis_dict); print(len_Arabidopsis_degs); print(Mod_Size_Percent);

outfile2 = open('Arabidopsis_DEGMod_Enrichment.txt','w')
outfile2.write('Module'+'\t'+'Observed_DEGs'+'\t'+'Expected_DEGs'+'\t'+'Fold_Change'+'\t'+'Fisher_Test'+'\n')
Arabidopsis_degs_count = {}; Arabidopsis_degs_counter = []
for j in range(Total_Modules):
	Arabidopsis_count = sum(x == j for x in Arabidopsis_dict.values()); #print(Arabidopsis_count);
	Arabidopsis_degs_count[j] = Arabidopsis_count	# dictionary where key = modules & values = DEG count for Pathogenesis
	Arabidopsis_degs_counter.append(Arabidopsis_count)	# List of DEG counts in each module for Pathogenesis
	Arabidopsis_deg_mod = float(Arabidopsis_count)	# DEG in each Module
	Arabidopsis_nondeg_mod = float(Mod_Size[j])-float(Arabidopsis_count)	# Non-DEG in each Module
	Arabidopsis_deg_nonmod = float(len_Arabidopsis_degs)-float(Arabidopsis_count)	# DEG not in each Module
	Arabidopsis_nondeg_nonmod = float(len_nodes)-float(Arabidopsis_nondeg_mod)	# Non-DEG not in each Module
	oddsratio, pvalue = stats.fisher_exact([[Arabidopsis_deg_mod, Arabidopsis_nondeg_mod], [Arabidopsis_deg_nonmod, Arabidopsis_nondeg_nonmod]]); # p-values of DEGs in each Module
	Arabidopsis_observed_degs_mod = Arabidopsis_count;	# DEGs observed in each Module
	Arabidopsis_expected_degs_mod = (len_Arabidopsis_degs*Mod_Size_Percent[j])/100	# DEGs in each Module, expected by random chance
	Arabidopsis_fold_change = float(Arabidopsis_observed_degs_mod)/float(Arabidopsis_expected_degs_mod)	# Fold Change of DEGs in each Module
	#print("P-value = ", pvalue);
	outfile2.write(str(j)+'\t'+str(Arabidopsis_observed_degs_mod)+'\t'+str(Arabidopsis_expected_degs_mod)+'\t'+str(Arabidopsis_fold_change)+'\t'+str(pvalue)+'\n')



