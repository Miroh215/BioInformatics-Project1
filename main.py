import sys

amino_dict = {
    "M": ["AUG"],  # Methionine (Start Codon)
    "A": ["GCU", "GCC", "GCA", "GCG"],  # Alanine
    "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],  # Arginine
    "N": ["AAU", "AAC"],  # Asparagine
    "D": ["GAU", "GAC"],  # Aspartic acid (Aspartate)
    "C": ["UGU", "UGC"],  # Cysteine
    "Q": ["CAA", "CAG"],  # Glutamine
    "E": ["GAA", "GAG"],  # Glutamic acid (Glutamate)
    "G": ["GGU", "GGC", "GGA", "GGG"],  # Glycine
    "H": ["CAU", "CAC"],  # Histidine
    "I": ["AUU", "AUC", "AUA"],  # Isoleucine
    "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],  # Leucine
    "K": ["AAA", "AAG"],  # Lysine
    "F": ["UUU", "UUC"],  # Phenylalanine
    "P": ["CCU", "CCC", "CCA", "CCG"],  # Proline
    "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],  # Serine
    "T": ["ACU", "ACC", "ACA", "ACG"],  # Threonine
    "W": ["UGG"],  # Tryptophan
    "Y": ["UAU", "UAC"],  # Tyrosine
    "V": ["GUU", "GUC", "GUA", "GUG"],  # Valine
    "B": ["GAU", "GAC"],  # Aspartic acid or Asparagine
    "X": ["UAA", "UAG", "UGA"],  # Stop Codon
}


def translation(strand: str):
    if 'u' in strand:
        pass


def get_file():
    in1 = input("Input the filename for the first sequence: ")
    try:
        with open(in1, 'r') as f:
            pass
    except FileNotFoundError:
        print("Invalid file.")
        exit()

    in2 = input("Input the filename for the second sequence: ")
    try:
        with open(in2, 'r') as f:
            pass
    except FileNotFoundError:
        print("Invalid file.")
        exit()

    return in1, in2


def extract_string(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        out_string = lines[1]
    return out_string


def extract_matrix(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        letters = lines[1].strip().split(",")
        mat = {}
        for letter in letters:
            mat[letter] = []
        for i in range(2, len(lines[2:])+2):
            temp = lines[i].strip().split(",")
            for x in range(len(temp)):
                mat[letters[x]].append(int(temp[x]))
    return mat, letters


def get_score(mat: dict, o: list, combo: tuple):
    return mat[combo[0]][o.index(combo[1])]


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        s1 = extract_string("Sample-Inputs/sequenceA1.txt")
        s2 = extract_string("Sample-Inputs/sequenceA2.txt")
        matrix, order = extract_matrix("Files/BLOSUM62.txt")
        # print(get_score(matrix, order, ("A", "C")))
        # print(get_score(matrix, order, ("C", "A")))
        # print(get_score(matrix, order, ("Y", "E")))
    else:
        infile, infile2 = get_file()
        s1 = extract_string(infile)
        s2 = extract_string(infile2)

    print(s1)
    print(s2)
