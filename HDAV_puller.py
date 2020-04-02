#Library
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from collections import namedtuple
from operator import methodcaller
from collections import defaultdict
import glob

#DICTIONARIES
genes_HDAVs = defaultdict(set)	#Number of HDAVs per gene across all species across all variants
spp_HDAVs = defaultdict(set) #Number of HDAVs per species across all variants found in all species
found_HDAVs = set() #HDAVs found in any species
spp_var_gaps = defaultdict(lambda: defaultdict(int))
spp_var_HDAVs = defaultdict(lambda: defaultdict(int))
var_info = defaultdict(str) # line with input data about variant. This will be used for output

#INPUT: Pathogen file containing HDAVs
varFile=open("output_LikelyPathogenic.txt", 'r')
#varFile=open("output_Pathogenic.txt", 'r')
speciesListFile = open("species_list.txt",'r')

#Make list of species, so species is always outputed in the same order
species_list = list()

for speciesLine in speciesListFile:
	species_list.append(speciesLine.rstrip())

for varLine in varFile:
	codonPos = 0	#position in codon that the variant is in (0,1,2). Set to -1 until it is determined
	varSplit=varLine.rstrip().split()
	varId = varSplit[0]
	if varSplit[10] == "None":
		continue	
	geneId = varSplit[2]
	if geneId == "Gene_ID":
		continue		
	#print(geneId)

	var_info[varId] = varLine.rstrip()
	refCodon = varSplit[8]
	#print(refCodon)
	
	#figuring out codon position
	for l in refCodon:
		if l.isupper():
			break	
		else:
			codonPos+=1
	#print(codonPos)
	
	#translating to find ref AA
	refAA=Seq(refCodon,generic_dna).translate()
	#print(refAA)

	#translating to find variant AA	
	varCodon = varSplit[9]
	varAA = Seq(varCodon,generic_dna).translate()
	#print(varAA)

	alnPos = int(varSplit[10])
	
	filegrabber = "/home/rumika.mascarenhas/PRANK_OUT/Processed_Prank/*"+geneId+"*"
	for filename in glob.glob(filegrabber):
		#print(filename)
		sppFile = open(filename,'r')
		for sppLine in sppFile:
			sppLine = sppLine.rstrip()
			if ">" in sppLine:
				species = sppLine[1:]
			else:
				#print(sppLine)
				
				#Pull out codon for specific species
				if codonPos == 0:
					sppCodon = sppLine[alnPos-1:alnPos+2]
				elif codonPos == 1:
					sppCodon = sppLine[alnPos-2:alnPos+1]
				elif codonPos == 2:
					sppCodon = sppLine[alnPos-3:alnPos]
				#print(sppCodon)
				if '-' in sppCodon:
						spp_var_gaps[varId][species] += 1
						continue	

				sppAA = Seq(sppCodon,generic_dna).translate()
				#print(sppAA)
			
				#check to see if HDAV is present in species
				if varAA == sppAA:
					#print(varId,geneId,species,varAA,sppAA)
					genes_HDAVs[geneId].add(varId)
					spp_HDAVs[species].add(varId)
					found_HDAVs.add(varId)
					spp_var_HDAVs[varId][species] += 1
		sppFile.close()
	#break	
HDAVs_per_gene = open("HDAVs_per_gene_LP.txt", "w")
for geneKey in genes_HDAVs:
	print(geneKey,len(genes_HDAVs[geneKey]))
	HDAVs_per_gene.write(geneKey+"\t"+str(len(genes_HDAVs[geneKey]))+"\n")
HDAVs_per_gene.close()

HDAVs_per_species = open("HDAVs_per_species_LP.txt","w")
for sppKey in spp_HDAVs:
	print(sppKey,len(spp_HDAVs[sppKey]))
	HDAVs_per_species.write(sppKey+"\t"+str(len(spp_HDAVs[sppKey]))+"\n")
HDAVs_per_species.close()

print(len(found_HDAVs))

HDAVs_found_file = open("HDAVs_in_nonhumans_LP.txt","w")
for HDAV in found_HDAVs:
	HDAVs_found_file.write(HDAV + "\n")
HDAVs_found_file.close()

#print(found_HDAVs)
##master spreadsheet
master_spreadsheet = open("likely_pathogenic_master_spreadsheet.txt","w")
master_spreadsheet.write("varID\tLocation\tGene_ID\tTranscript_ID\tRef_Var\tAlt_Var\tRef_Codon\tAlt_Codon\tRef_AA\tAlt_AA\tAlignment_Position\t")
for spp in species_list:
	master_spreadsheet.write(spp+"\t")
master_spreadsheet.write("\n")

for varKey in var_info:
	master_spreadsheet.write(var_info[varKey]+"\t")
	for spp in species_list:
		master_spreadsheet.write(str(spp_var_HDAVs[varKey][spp])+"\t")
	master_spreadsheet.write("\n")
master_spreadsheet.close()
