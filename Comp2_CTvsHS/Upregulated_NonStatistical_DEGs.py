import re

file1 = open('Arabidopsis_DEGMod_Enrichment.txt','r')
for i in range(2):
	line1 = file1.readline()
module, observed, expected, foldchange, fishertest = [],[],[],[],[]
stats_module, stats_foldchange, stats_fishertest = [],[],[]
while line1:
	split_line1 = (line1.rstrip()).split('\t')
	if (float(split_line1[3]) > 1): # DEG Enriched Modules with Fold Change is greater than 1
		module.append(split_line1[0])
		observed.append(split_line1[1])
		expected.append(split_line1[2])
		foldchange.append(split_line1[3])
		fishertest.append(split_line1[4])
		if (float(split_line1[4]) <= 0.05): # DEG Enriched Modules with p-value is less than equal to 0.05
			stats_module.append(split_line1[0])
			stats_foldchange.append(split_line1[3])
			stats_fishertest.append(split_line1[4])
	line1 = file1.readline()

file5 = open('WGCNA_Arabidopsis_Modules_Clusters.txt','r')
m=0; mod, transcript = [],[];
for line5 in file5:
	split_line5 = (line5.rstrip()).split('\t')
	mod.append(str(m))
	transcript.append(split_line5)
	m += 1

file2 = open('Comp2_CTvsHS.txt_downregulated.txt','r')
line2 = file2.readline()
locusid_down,transcriptid_down,genename_down,foldchanges_down,pvalue_down,koid_down,cogid_down,cogcategory_down = [],[],[],[],[],[],[],[];
while line2:
	split_line2 = (line2.rstrip()).split('\t')
	locusid_down.append(split_line2[0])
	#transcriptid_down.append(split_line2[1])
	#genename_down.append(split_line2[2])
	foldchanges_down.append(split_line2[3])
	#pvalue_down.append(split_line2[5])
	#koid_down.append(split_line2[5])
	#cogid_down.append(split_line2[6])
	#cogcategory_down.append(split_line2[7])
	line2 = file2.readline()

file3 = open('Comp2_CTvsHS.txt_upregulated.txt','r')
for i in range(2):
	line3 = file3.readline()
locusid, transcriptid, genename, foldchanges, pvalue, koid, cogid, cogcategory = [],[],[],[],[],[],[],[];
while line3:
	split_line3 = (line3.rstrip()).split('\t')
	locusid.append(split_line3[0])
	#transcriptid.append(split_line3[1])
	#genename.append(split_line3[2])
	foldchanges.append(split_line3[3])
	#pvalue.append(split_line3[5])
	#koid.append(split_line3[5])
	#cogid.append(split_line3[6])
	#cogcategory.append(split_line3[7])
	line3 = file3.readline()

file4 = open('Comp2_CTvsHS.txt','r')
for i in range(2):
	line4 = file4.readline()
locus_id, gene_name, fold_change, p_value, ko_id, cog_id, cog_category = [],[],[],[],[],[],[];
while line4:
	split_line4 = (line4.rstrip()).split('\t')
	locus_id.append(split_line4[0])
	#gene_name.append(split_line4[1])
	fold_change.append(split_line4[2])
	p_value.append(split_line4[6])
	#ko_id.append(split_line4[6])
	#cog_id.append(split_line4[7])
	#cog_category.append(split_line4[8])
	line4 = file4.readline()
"""
file6 = open('Arabidopsis_Annotations.tsv','r')
for i in range(2):
	line6 = file6.readline()
locus_tag, transcript_tag = [],[];
while line6:
	split_line6 = (line6.rstrip()).split('\t')
	locus_tag.append(split_line6[6])
	transcript_tag.append(split_line6[7]); #print(split_line6[6], split_line6[7]);
	line6 = file6.readline()
"""

outfile1 = open('All_Transcripts_In_DEGenrichedModules.txt','w')
outfile1.write('Module'+'\t'+'Transcript_Id'+'\t'+'FoldChange_Module'+'\t'+'FisherTest_Module'+'\n')
outfile2 = open('UpRegulatedNonDEGs_In_DEGenrichedModules.txt','w')
outfile2.write('Module'+'\t'+'Locus_Tag'+'\t'+'Expression_Fold_Change'+'\t'+'Expression_Pvalue'+'\t'+'FoldChange_Module'+'\t'+'FisherTest_Module'+'\n')
outfile3 = open('DownRegulatedNonDEGs_In_DEGenrichedModules.txt','w')
outfile3.write('Module'+'\t'+'Locus_Tag'+'\t'+'Expression_Fold_Change'+'\t'+'Expression_Pvalue'+'\t'+'FoldChange_Module'+'\t'+'FisherTest_Module'+'\n')
modid, transids = [],[];
for x in range(len(module)): # DEG Enriched Modules
	if str(module[x]) in mod:
		index1 = int(mod.index(str(module[x])))
		mod_transcripts = transcript[index1]; 
		for t in range(len(mod_transcripts)): # All transcripts in DEG Enriched Modules
			outfile1.write(str(module[x])+'\t'+mod_transcripts[t]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\n')
			modid.append(str(module[x]))
			transids.append(mod_transcripts[t])
			#if (mod_transcripts[t] in transcript_tag):
			#	index2 = int(transcript_tag.index(mod_transcripts[t]))
			locus_gene = mod_transcripts[t]
			if (locus_gene in locus_id):
				index3 = int(locus_id.index(locus_gene))
				if ((fold_change[index3]) != 'NA'):
					if (float(fold_change[index3])>1): # Up-regulated DEGs. You may need to change to >=4 instead of >2 
						if (locus_gene not in locusid): # Up-regulated DEGs not in paper.
							outfile2.write(str(module[x])+'\t'+locus_gene+'\t'+fold_change[index3]+'\t'+p_value[index3]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\n')
					if (float(fold_change[index3])<-1): # Down-regulated DEGs 
						if (locus_gene not in locusid_down): # Down-regulated DEGs not in paper.
							outfile3.write(str(module[x])+'\t'+locus_gene+'\t'+fold_change[index3]+'\t'+p_value[index3]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\n')











