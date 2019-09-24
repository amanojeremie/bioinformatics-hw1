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
@param regionLen length of the region to look at
"""
def getHdyroRegion(aminoAcid, regionLen):
    acid = aminoAcid.split('-') #remove "-"
    acid.remove('STOP')
    results = [] #store potential region
    for i in range(0, len(acid) - regionLen + 1): #loop through amino acid with given length
        hydroSum = 0
        region = acid[i : i + regionLen] #get subregion 
        for j in range(0, regionLen): #sum up the hydro. 
            hydroSum += getHydro(region[j])
        if hydroSum >= 25: # if the sum is greater than a critical value, then it is a transport potein
            results.append(region) # store subregion
    print(results, end="\n\n") #print subregion. Will decide how to highlight the region

"""
get the individual hydrophobicity
@param acid individual aminio acid
@return  kdHydrophobicity by default
"""
def getHydro(acid): #get individual hydro. given an acid
    for acidHydro in hydrophobicity:
        if acid in acidHydro[1]:
            
            return acidHydro[2]



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
            
			getHdyroRegion(str(translation), 17)
           
			
if __name__ == "__main__":
	main(sys.argv)