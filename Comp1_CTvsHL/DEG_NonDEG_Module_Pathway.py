import re

file1 = open('WGCNA_Arabidopsis_Modules_Clusters.txt','r')
mod_count = 0; module, transcripts, module_length = [],[],[];
for line1 in file1:
	split_line1 = (line1.rstrip()).split('\t')
	module.append(mod_count)
	transcripts.append(split_line1)
	module_length.append(len(split_line1))
	mod_count += 1

file2 = open('Comp1_CTvsHL_DEGs.txt','r')
for i in range(2):
	line2 = file2.readline()
transcript_id, locus_tag, fold_change, gene_name, ko_id, pvalue, cog = [],[],[],[],[],[],[];
while line2:
	split_line2 = (line2.rstrip()).split('\t')
	transcript_id.append(split_line2[0])
	fold_change.append(split_line2[3])
	line2 = file2.readline()

file3 = open('UpRegulatedNonDEGs_In_DEGenrichedModules.txt','r')
for i in range(2):
	line3 = file3.readline()
modules, transcriptid, locustag, genename, expression_foldchange, expression_pvalue = [],[],[],[],[],[];
koid, cogid, module_foldchange, module_fishertest, cog_category = [],[],[],[],[];
while line3:
	split_line3 = (line3.rstrip()).split('\t')
	modules.append(split_line3[0])
	transcriptid.append(split_line3[1]); #print(transcriptid);
	expression_foldchange.append(split_line3[2])
	expression_pvalue.append(split_line3[3])
	module_foldchange.append(split_line3[4])
	module_fishertest.append(split_line3[5])
	line3 = file3.readline()

file4 = open('DownRegulatedNonDEGs_In_DEGenrichedModules.txt','r')
for i in range(2):
	line4 = file4.readline()
while line4:
	split_line4 = (line4.rstrip()).split('\t')
	modules.append(split_line4[0])
	transcriptid.append(split_line4[1]); #print(transcriptid);
	expression_foldchange.append(split_line4[2])
	expression_pvalue.append(split_line4[3])
	module_foldchange.append(split_line4[4])
	module_fishertest.append(split_line4[5])
	line4 = file4.readline()

outfile = open('DEG_NonDEG_Module_Pathway.txt','w')
outfile.write('Module'+'\t'+'DEG_Count'+'\t'+'Up&DownRegulated_NonDEG_Count'+'\t'+'Total_DEG_Count'+'\t'+'Module_Gene_Count'+'\t'+'DEG_Percent_In_Module'+'\t'+'DEG_Pval_Percent_In_Module'+'\n')
for i in range(len(transcripts)):
	transcript_deg_count = 0; transcript_nondeg_count = 0;
	for j in range(len(transcripts[i])):
		if ((transcripts[i])[j] in transcript_id):
			transcript_deg_count += 1
		if ((transcripts[i])[j] in transcriptid):
			transcript_nondeg_count += 1
	transcript_count = transcript_deg_count+transcript_nondeg_count
	#print(module[i], transcript_deg_count, transcript_nondeg_count, transcript_count, module_length[i])
	#print(module[i], transcript_deg_count, transcript_nondeg_count, transcript_count, module_length[i], str(float(transcript_count/int(module_length[i]))))
	outfile.write(str(module[i])+'\t'+str(transcript_deg_count)+'\t'+str(transcript_nondeg_count)+'\t'+str(transcript_count)+'\t'+str(module_length[i])+'\t'+str(float(transcript_count/int(module_length[i])))+'\t'+str(float(transcript_deg_count/int(module_length[i])))+'\n') # DEG_FC4_Pval_Percent_In_Module: These are the DEGs that show Fold change >= 4 and, Pvalue < 0.01.




