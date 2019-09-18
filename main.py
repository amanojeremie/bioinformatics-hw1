import sys
from translate import translate

def main(argv):
	"""
		Given a FASTA file, translate from RNA to amino acids
		
		Usage:
		python main.py fasta.txt
	"""
	if(len(argv) < 2):
		print("Not enough arguments")
		return
	
	file = open(argv[1], "r")
	for line in file:
		if line[0] == ">":
			name = line[1:-1] #Omits the > and the newline character
			sequence = file.readline()[:-1] #Omits the newline character
			print(name)
			print(translate(sequence))
if __name__ == "__main__":
	main(sys.argv)