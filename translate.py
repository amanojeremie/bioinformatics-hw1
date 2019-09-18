"""
	Amino acids are represented as a 4-tuple
	The first element is the full name of the acid
	The second element is the short representation of the acid
	The third element is the letter representation of the acid
	The fourth element is a list of accepted RNA translations for the acid
"""
AMINO_ACIDS = [
	("Alanine", "Ala", "A", ["GCU", "GCC", "GCA", "GCG"]),
	("Arginine", "Arg", "R", ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"]),
	("Asparagine", "Asn", "N", ["AAU", "AUC"]),
	("Aspartic Acid", "Asp", "D", ["GAU", "GAC"]),
	("Cysteine", "Cys", "C", ["UGU", "UGC"]),
	("Glutamic Acid", "Glu", "E", ["GLU", "GAG"]),
	("Glutamine", "Gln", "Q", ["CAA", "CAG"]),
	("Glycine", "Gly", "G", ["GGU", "GGC", "GGA", "GGG"]),
	("Histidine", "His", "H", ["CAU", "CAC"]),
	("Isoleucine", "Ile", "I", ["AUU", "AUC", "AUA"]),
	("Leucine", "Leu", "L", ["CUU", "CUC", "CUA", "CUG"]),
	("Lysine", "Lys", "K", ["AAA", "AAG"]),
	("Methionine", "Met", "M", ["AUG"]),
	("Phenylalanine", "Phe", "F", ["UUU", "UUC"]),
	("Proline", "Pro", "P", ["CCU", "CCC", "CCA", "CCG"]),
	("Serine", "Ser", "S", ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"]),
	("Threonine", "Thr", "T", ["ACU", "ACC", "ACA", "ACG"]),
	("Tyrosine", "Tyr", "Y", ["UAU", "UAC"]),
	("Valine", "Val", "V", ["GUU", "GUC", "GUA", "GUG"]),
	("STOP", "STOP", "STOP", ["UAA", "UAG", "UGA"])
]

class AcidSequence:
	"""
		Class that represents a translated Acid sequence given from
		an RNA sequence.
	"""
	def __init__(self, acidArray):
		self.acidArray = acidArray

	def __str__(self):
		"""
			Creates a string representation of the sequence using the short representation
		"""
		str = ""
		for aminoAcid in self.acidArray:
			str += aminoAcid[1] + "-"
		return str[:-1] #Omits the last dash

def tripletToAcid(triplet):
	"""
		Converts a triplet of characters into an amino acid
	"""
	for aminoAcid in AMINO_ACIDS:
		if triplet in aminoAcid[3]:
			return aminoAcid

def translate(str):
	"""
		Translates an RNA sequence string into an amino acid sequence.
		Returns an AcidSequence
	"""

	if len(str) % 3 != 0: #Sequence length has to be divisible by 3
		print("Not translatable")
		return
	translation = []
	for i in range(0, len(str) // 3): #From 0 to number of triplets
		triplet = str[i * 3 : i * 3 + 3] #Extracts a triplet from the string
		acid = tripletToAcid(triplet)
		if acid != None:
			if acid[0] == "STOP": #If a stop codon is reached, stop translating
				continue
			translation.append(acid)
	return AcidSequence(translation)
