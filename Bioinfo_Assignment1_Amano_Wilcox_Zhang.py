"""
	Amino acids are represented as a 4-tuple
	The first element is the full name of the acid
	The second element is the 3-letter representation of the acid
	The third element is the 1-letter representation of the acid
	The fourth element is a list of accepted RNA codons for the acid
"""
AMINO_ACIDS = [
	("Alanine", "Ala", "A", ["GCU", "GCC", "GCA", "GCG"]),
	("Arginine", "Arg", "R", ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"]),
	("Asparagine", "Asn", "N", ["AAU", "AAC"]),
	("Aspartic Acid", "Asp", "D", ["GAU", "GAC"]),
	("Cysteine", "Cys", "C", ["UGU", "UGC"]),
	("Glutamic Acid", "Glu", "E", ["GAA", "GAG"]),
	("Glutamine", "Gln", "Q", ["CAA", "CAG"]),
	("Glycine", "Gly", "G", ["GGU", "GGC", "GGA", "GGG"]),
	("Histidine", "His", "H", ["CAU", "CAC"]),
	("Isoleucine", "Ile", "I", ["AUU", "AUC", "AUA"]),
	("Leucine", "Leu", "L", ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"]),
	("Lysine", "Lys", "K", ["AAA", "AAG"]),
	("Methionine", "Met", "M", ["AUG"]),
	("Phenylalanine", "Phe", "F", ["UUU", "UUC"]),
	("Proline", "Pro", "P", ["CCU", "CCC", "CCA", "CCG"]),
	("Serine", "Ser", "S", ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"]),
	("Threonine", "Thr", "T", ["ACU", "ACC", "ACA", "ACG"]),
	("Tryptophan", "Trp", "W", ["UGG"]),
	("Tyrosine", "Tyr", "Y", ["UAU", "UAC"]),
	("Valine", "Val", "V", ["GUU", "GUC", "GUA", "GUG"]),
	("STOP", "STOP", "STOP", ["UAA", "UAG", "UGA"])
]

class AminoAcidSequence:
	"""
		Class that represents a translated Amino Acid sequence given from
		an RNA sequence.
	"""
	def __init__(self, acidArray):
		self.acidArray = acidArray

	def __str__(self):
		"""
			Creates a string representation of the amino acid sequence using the 1-letter representation
		"""
		str = ""
		for aminoAcid in self.acidArray:
			str += aminoAcid[2]
		return str

def codonToAminoAcid(codon):
	"""
		Converts a codon into an amino acid
	"""
	for aminoAcid in AMINO_ACIDS:
		if codon in aminoAcid[3]: #If codon is in the list of accepted codons
			return aminoAcid
	print("UNKNOWN:", codon) #If no match found, print an unkwown warning

def translate(str):
	"""
		Translates an RNA sequence string into an amino acid sequence.
		Returns an AminoAcidSequence
	"""

	if len(str) % 3 != 0: #Sequence length has to be divisible by 3
		print("Not translatable")
		return

	translation = []
	for i in range(0, len(str) // 3): #From 0 to number of codons
		codon = str[i * 3 : i * 3 + 3] #Extracts a codon from the string
		acid = codonToAminoAcid(codon)
		if acid != None:
			translation.append(acid)
			if acid[0] == "STOP": #If a stop codon is reached, stop translating
				break
	return AminoAcidSequence(translation)

import os, sys
from translate import translate

if sys.platform.lower() == "win32": #Enables color output in Windows CLIs, not needed for other OSs
	os.system('color')

"""
First element is three letter representation of amino acid
Second element: one letter representaion of amino acid
Third element: kdHydrophobicity
Fourth element: lHydro
"""

hydrophobicity = [
	("Ile", "I", 4.5, 1),
	("Val", "V", 4.2,  1),
	("Leu", "L", 3.8, 1),
	("Phe", "F", 2.8, 1),
	("Cys", "C", 2.5, 1),
	("Met", "M", 1.9, 1),
	("Ala", "A", 1.8, 0),
	("Gly", "G", -0.4, -1),
	("Thr", "T", -0.7, -1),
	("Ser", "S", -0.8, -1),
	("Trp", "W", -0.9, 1),
	("Tyr", "Y", -1.3, 0),
	("Pro", "P", -1.6, -1),
	("His", "H", -3.2, 0),
	("Glu", "E", -3.5, -1),
	("Gln", "Q", -3.5, -1),
	("Asp", "D", -3.5, -1),
	("Asn", "N", -3.5, -1),
	("Lys", "K", -3.9, -1),
	("Arg", "R", -4.5, -1)
]

"""
get the total hydrophobicity given the region length
print out the regions that have high hydrophobicity
@param aminoAcid 
@param minLen minimal length of region
@param maxLen maximum length of region
@param cutoff the cutoff hydrophobicity to determine if a region is a hydrophobic
@param useLHydro A boolean that specifies if lHydro will be used instead of kdHydro
"""
def getHydrophobicRegions(aminoAcid, minLen, maxLen, cutoff, useLHydro=False):
	acid = aminoAcid[0 : -4] # remove "stop" codon
	results = [] #store potential region
	resultsIndex = [] #store region index
	preI = -999
	preLen = 0
	for i in range(0, len(acid) - minLen + 1): #loop through amino acid with given length
		for j in range(minLen, maxLen): #try different length of region: from min to max
			hydroSum = 0  
			if((i + j) < len(acid)): #avoid out of bound
				region = acid[i : i + j] #get subregion     
				for k in range(0, j ): #sum up the hydrophobicity values
					hydroSum += getHydro(region[k], useLHydro)
					flag = (i - preI) >= preLen #the start index of second region should at least j away from the start of previous region
					if (hydroSum >= cutoff) and flag: # if the sum is greater than a critical value, then it is a transport potein  
						
						charge = compareCharges(i, i + j, acid, 30) #Compares which side of the region has more postively charged amino acids
						if charge == -1: 
							region = "+" + region
						elif charge == 1:
							region = region + "+"

						results.append(region) # store subregion
						preLen = len(region)
						preI = i
						index = (i, i + j)
						resultsIndex.append(index)
						break #once a region with high hydrophobicity is found stop increasing the length of region                    
					else:
						continue #continue increase the size of region

	#Creates a string where the transmembrane domains in the amino acid sequence are highlighted
	coloredAminoAcidSequence = ""
	lastIndex = 0
	for index in resultsIndex:
		coloredAminoAcidSequence = coloredAminoAcidSequence + acid[lastIndex : index[0]] + "\x1b[6;30;42m"
		coloredAminoAcidSequence = coloredAminoAcidSequence + acid[index[0] : index[1]] + "\x1b[0m"
		lastIndex = index[1]
	coloredAminoAcidSequence = coloredAminoAcidSequence + acid[lastIndex:]

	print("Amino acid sequence:", coloredAminoAcidSequence)   
	print("Transmembrane domain sequence indices:", resultsIndex, end="\n")
	print("Transmembrane domain sequences:", results, end="\n\n") #print subregion.

"""
get the individual hydrophobicity
@param acid individual amino acid
@return  hhHydrophobicity by default
"""
def getHydro(acid, useLHydro): #get individual hydrophobicity given an acid
	for acidHydro in hydrophobicity:
		if acid in acidHydro[1]:            
			return acidHydro[3] if useLHydro else acidHydro[2]

"""
Compares the charges of two sides of a region in an amino acid sequence to determine which side has more charged amino acids
@param start the start index of the region
@param end the end index of the region
@param acidSeq the amino acid sequence
@param size the window size to check on both sides
@return 0 if no side, -1 if more charges before the region, 1 if more charges after the region
"""
#Array of positively charged amino acids
charged = ["R", "H", "K"]
def compareCharges(start, end, acidSeq, size):
	if start - size < 0 or end + size >= len(acidSeq): #Insure within bounds
		return 0
	
	leftChargeSum = 0
	for acid in acidSeq[start - size : start]:
		if acid in charged: #If amino acid is positively charged, increment sum
			leftChargeSum += 1
	
	rightChargeSum = 0
	for acid in acidSeq[end : end + size]:
		if acid in charged:
			rightChargeSum += 1
	if leftChargeSum < rightChargeSum:
		return 1
	elif rightChargeSum < leftChargeSum:
		return -1
	else:
		return 0

def main(argv):
	"""
		Given a FASTA file, translate from RNA to amino acids
		
		Usage:
		python main.py fasta.txt
	"""
	if(len(argv) < 2):
		print("Given an RNA sequence, predicts possible transmembrane domains within the sequence")
		print("Usage: python", argv[0], "fasta.txt (-l)")
		print("Options:")
		print("\t-l:", "Use simplified hydrophobicity table")
		return
	
	file = open(argv[1], "r")
	for line in file:
		if line[0] == ">":
			name = line[1:-1] #Omits the > and the newline character
			sequence = file.readline()[:-1] #Omits the newline character
			print(name)
			translation = translate(sequence)
			
			if len(argv) >= 3 and argv[2] == "-l":
				getHydrophobicRegions(str(translation), 20, 25, 6, True)
			else:
				getHydrophobicRegions(str(translation), 20, 25, 30)

		   
			
if __name__ == "__main__":
	main(sys.argv)