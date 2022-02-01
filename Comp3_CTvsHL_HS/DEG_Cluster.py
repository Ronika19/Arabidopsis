file1 = open('Comp3_CTvsHL_HS.txt','r')
file4 = open('Arabidopsis_Cluster_DEG_Comp3_CTvsHL_HS.txt','w')

for i in range(2):
	line1 = file1.readline()
while line1:
	line1 = line1.rstrip()
	split_line1 = line1.split('\t')
	genes1 = split_line1[0]
	fold_change1 = split_line1[2]
	#if ((float(fold_change1) >= 1) and (float(split_line1[6]) < 0.05)):
	if (split_line1[6] != 'NA'):
		if (float(split_line1[6]) < 0.05):
			file4.write(genes1+"\t")
	#if ((float(fold_change1) <= -1) and (float(split_line1[6]) < 0.05)):
	#	file4.write(genes1+"\t")
	line1 = file1.readline()
file4.write("\n")
file4.close()



