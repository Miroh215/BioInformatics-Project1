

amino_dict = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic acid (Aspartate)",
    "C": "Cysteine",
    "Q": "Glutamine",
    "E": "Glutamic acid (Glutamate)",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
    "B": "Aspartic acid or Asparagine",
    "Z": "Glutamine or Glutamic acid.",
    "X": "Any amino acid.",
    " ": "termination codon"
}


def get_form(failed=False):
    if failed:
        print("Incorrect format. Please try again.")
    else:
        print("Will the sequence files be for nucleotides or peptides?")
    return input("Input 'n' for nucleotides, 'p' for peptides, or 'q' to quit: ").strip()[0]


def get_filename(failed=False):
    if failed:
        print("Please ensure the input files are in proper FASTA format and/or exist and try again.")
    return input("Input the filename: ")


def translation(strand: str):
    if 'u' in strand:
        pass


if __name__ == "__main__":
    form = get_form()
    while form not in ['p', 'n', 'q']:
        form = get_form(failed=True)

    if form == 'q':
        exit()
    elif form == 'n':
        print("You have selected nucleotide.")
    else:
        print("You have selected peptide.")

    invalid = True

    while invalid:
        infile = get_filename()
        if infile == 'q':
            exit()
        try:
            f = open(infile, 'r')
            invalid = False
        except FileNotFoundError:
            print("Invalid file.")

    f.close()
    print("This file isn't finished yet.")
