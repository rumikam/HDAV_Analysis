#get the CDS from the ensemble output file for each HDAV, and
#gets the position in the human transcript
#updated: March 19th, 2020
filename = "VEP_Pathogenic.txt"
ensEMBL_files = "aligned_files.txt"
new_pos = {}

output = open("output_Pathogenic.txt", "w")
output_error = open("output_error_Pathogenic.txt", "w")
alignment_length_error = open("output_alignment_length_errorP.txt", "w")
NotInEnsEMBL = open("output_NotInEnsEMBL_P.txt", "w")
NotCanonicalSNV = open("output_NotCanonicalSNV_P.txt", "w")


def isCanonicalSNV(canonical, CDS_pos):
	if((canonical == "YES") and ("-" not in CDS_pos)):
		return True
	else:
		return False

def isInEnsEMBL(Gene):
	with open(ensEMBL_files, 'r') as y: #EnsEMBL_files is a txt with all the EnsEMBL files
			for ensEMBL_geneid in y:
				if Gene in ensEMBL_geneid:
					return ensEMBL_geneid

def getNewPos(file, CDS_pos):
	with open("/home/rumika.mascarenhas/PRANK_OUT/Processed_Prank/" + file.strip("\n"), 'r') as z:
		next(z)
		for fasta in z:
		        #print(fasta)
			ClinVar_CDS_pos = int(CDS_pos)
			count = 0 #increments when a gap occurs
			if("Homo_sapiens" in fasta):
				align = str(next(z))
				align = list(align)
                                print(align)
				if(len(align) >= ClinVar_CDS_pos):
					for i in range(0,(ClinVar_CDS_pos-1)):
						if(align[i] == "-" ):
							count += 1
				else:
					alignment_length_error.write(location + "\t" + gene + "\t" + transcript + "\n")
			new_position = str(ClinVar_CDS_pos + count)
			print(count)
			print(new_position)
			return new_position

with open(filename, 'r') as x:
	next(x)
	for line in x:
	    parts = line.split()
	    ID = parts[0]
            #print(ID)
	    location = parts[1]
	    alt = parts[2]
	    #print(alt)
	    gene = parts[3]
	    transcript = parts[4]
	    ClinVar_CDS_pos = parts[8] #Coding sequence position
	    #Ref = parts[13]
            codon = parts[11]
            codon_split = codon.split("/")
	    ref_codon = codon_split[0]
            if(len(codon_split) == 2):
		    alt_codon = codon_split[1]
            else:
		    alt_codon = "-"	
	    for l in ref_codon:
	            if l.isupper():
			Ref = l
	    #print(Ref)
	    AA = parts[10]
	    #AA in format: Ref_AA/Alt_AA
	    AA_split = AA.split("/")
	    if(len(AA_split) == 2):
		    Ref_AA = AA_split[0]
		    Alt_AA = AA_split[1]

	    #symbol = parts[19]
	    if("CANONICAL=YES" in parts[13]):
		    Canonical = "YES"
	    else:
		    Canonical = "NO"

	    if(isCanonicalSNV(Canonical, ClinVar_CDS_pos)):
			file = isInEnsEMBL(gene)
			if(file != None):
				new_pos[ID + "\t" + location + "\t" + gene + "\t" + transcript + "\t" + Ref + "\t" + alt + "\t" + Ref_AA + "\t" + Alt_AA + "\t" + ref_codon + "\t" + alt_codon] = getNewPos(file, ClinVar_CDS_pos)
			if(file == None):
				NotInEnsEMBL.write(location + "\t" + gene + "\t" + transcript + "\n")
            else:
			NotCanonicalSNV.write(location + "\t" + gene + "\t" + transcript + "\n")

for j in new_pos.keys():
	output.write(j)
	output.write("\t")
	output.write(new_pos[j])
	output.write("\n")
