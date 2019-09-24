import sys
from translate import translate

"""
First element is three letter representation of amino acid
Second element: one letter representaion of amino acid
Third element: kdHydrophobicity
Fourth element: wwHydrophobicity
Fifth element: hhHydrophobicity
"""

hydrophobicity = [
    ("Ile", "I", 4.5, 0.31, -0.60 ),
    ("Val", "V", 4.2, -0.07, -0.31),
    ("Leu", "L", 3.8, 0.56, -0.55),
    ("Phe", "F", 2.8, 1.13, -0.32),
    ("Cys", "C", 2.5, 0.24, -0.13),
    ("Met", "M", 1.9, 0.23, -0.10),
    ("Ala", "A", 1.8, -0.17, 0.11),
    ("Gly", "G", -0.4, -0.01, 0.74),
    ("Thr", "T", -0.7, -0.14, 0.52),
    ("Ser", "S", -0.8, -0.13, 0.84),
    ("Trp", "W", -0.9, 1.85, 0.30),
    ("Tyr", "Y", -1.3, 0.94, 0.68),
    ("Pro", "P", -1.6, -0.45, 2.23),
    ("His", "H", -3.2, -0.96, 2.06),
    ("Glu", "E", -3.5, -2.02, 2.68),
    ("Gln", "Q", -3.5, -0.58, 2.36),
    ("Asp", "D", -3.5, -1.23, 3.49),
    ("Asn", "N", -3.5, -1.23, 3.49),
    ("Lys", "K", -3.9, -0.99, 2.71),
    ("Arg", "R", -4.5, -0.81, 2.58)
]

"""
get the total hydrophobicity given the region length
print out the regions that have high hydrophobicity
@param aminoAcid 
@param minLen minimal length of region
@param maxLen maximum length of region
"""
def getHdyroRegion(aminoAcid, minLen, maxLen):
    acid = aminoAcid[0 : len(aminoAcid) - 4] # remove "stop" codon
    print(acid)
    results = [] #store potential region
    resultsIndex = [] #store region index
    preI = -999
    preLen = 0
    for i in range(0, len(acid) - minLen + 1): #loop through amino acid with given length
        for j in range(minLen, maxLen): #try different length of region: from min to max
            hydroSum = 0  
            if((i + j) < len(acid)): #aviod out of bound
                region = acid[i : i + j] #get subregion     
                for k in range(0, j ): #sum up the hydro. 
                    hydroSum += getHydro(region[k])
                    flag = (i - preI) >= preLen #the start index of second region should at least j away from the start of previous region
                    if (hydroSum >= 45) and flag: # if the sum is greater than a critical value, then it is a transport potein  
                        results.append(region) # store subregion
                        preLen = len(region)
                        preI = i
                        index = (i, i + j)
                        resultsIndex.append(index)
                        break #once find a high hydro. stop increasing the length of region                    
                    else:
                        continue #continue increase the size of region
    print(resultsIndex, end="\n")
    print(results, end="\n\n") #print subregion. Will decide how to highlight the region

"""
get the individual hydrophobicity
@param acid individual aminio acid
@return  hhHydrophobicity by default
"""
def getHydro(acid): #get individual hydro. given an acid
    for acidHydro in hydrophobicity:
        if acid in acidHydro[1]:            
            return acidHydro[4]



def main(argv):
	"""
		Given a FASTA file, translate from RNA to amino acids
		
		Usage:
		python main.py fasta.txt
	"""
	if(len(argv) < 2):
		print("Translates RNA sequences in a FASTA file into an amino acid sequence")
		print("Usage: python", argv[0], "fasta.txt")
		return
	
	file = open(argv[1], "r")
	for line in file:
		if line[0] == ">":
			name = line[1:-1] #Omits the > and the newline character
			sequence = file.readline()[:-1] #Omits the newline character
			print(name)
			translation = translate(sequence)
            
			getHdyroRegion(str(translation), 17, 22)
           
			
if __name__ == "__main__":
	main(sys.argv)