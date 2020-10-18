"""
Author: Marcus Underwood
Class: CS490 Bioinformatics
Assignment: Project 1
Date: 10/14/2020 
"""


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
    if 'U' in strand:
        pairs = [strand[i:i+3] for i in range(0, len(strand), 3)]
    else:
        table = {
            "A": "A",
            "T": "U",
            "C": "C",
            "G": "G"
        }
        new_string = ""
        for char in strand:
            new_string += table[char]
        # print(f"Translation: {new_string}")
        pairs = [new_string[i:i+3] for i in range(0, len(new_string), 3)]
    output = ""
    for pair in pairs:
        for item in amino_dict:
            if pair in amino_dict[item]:
                if item == "X":
                    # print(pairs)
                    print(output)
                    return output
                else:
                    output += item
    # print(pairs)
    print(output)
    return output


def get_files():
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

    mat = input("Input the filename for the matrix: ")
    try:
        with open(mat, 'r') as f:
            pass
    except FileNotFoundError:
        print("Invalid file.")
        exit()

    return in1, in2, mat


def extract_string(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        out_string = lines[1].strip()
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


def semiglobal_align(string1, string2, mat, o, gap_pen):
    grid = [[None for x in range(len(string1) + 1)] for y in range(len(string2) + 1)]
    grid[0][0] = (0, 's')

    for i in range(1 ,len(string1)+1):
        grid[0][i] = (0, 'e')

    for i in range(1, len(string2)+1):
        grid[i][0] = (0, 'e')

    for x in range(len(string2)):
        for y in range(len(string1)):
            score = get_score(mat, o, (string1[y], string2[x]))
            best = max(grid[x][y][0] + score, grid[x+1][y][0] + gap_pen, grid[x][y+1][0] + gap_pen)
            if best == score + grid[x][y][0]:
                grid[x+1][y+1] = (best, 'd')
            elif best == grid[x+1][y][0] + gap_pen:
                grid[x+1][y+1] = (best, 'v')
            else:
                grid[x+1][y+1] = (best, 'h')

    for item in grid:
        out = ''
        for place in item:
            out += f'{place[0]:>8} '
        print(out)

    return grid


def semiglobal_backtrack(table, string1, string2):
    len1 = len(string1)
    len2 = len(string2)
    seq1 = ''
    seq2 = ''
    path = ''
    current = table[len2][len1]
    score = current[0]

    while current[1] != 's':
        if current[0] <= 0:
            break
        if current[1] == 'd':
            seq1 = string1[len1 - 1] + seq1
            seq2 = string2[len2 - 1] + seq2
            len1 -= 1
            len2 -= 1
            path = 'D' + path
        elif current[1] == 'h':
            seq1 = '-' + seq1
            seq2 = string2[len2 - 1] + seq2
            len2 -= 1
            path = 'H' + path
        elif current[1] == 'v':
            seq1 = string1[len1 - 1] + seq1
            seq2 = '-' + seq2
            len1 -= 1
            path = 'V' + path
        current = table[len2][len1]

    print(f"Alignment score: {score}")
    print(f"Forward path is: {path}")

    print("Aligned sequences:")
    print(seq1)
    print(seq2)


def global_align(string1, string2, mat, o, gap_pen):
    grid = [[None for x in range(len(string1)+1)] for y in range(len(string2)+1)]
    grid[0][0] = (0, 's')

    for i in range(1 ,len(string1)+1):
        grid[0][i] = (grid[0][i-1][0] + gap_pen, 'h')

    for i in range(1, len(string2)+1):
        grid[i][0] = (grid[i-1][0][0] + gap_pen, 'v')

    for x in range(len(string2)):
        for y in range(len(string1)):
            score = get_score(mat, o, (string1[y], string2[x]))
            best = max(grid[x][y][0] + score, grid[x+1][y][0] + gap_pen, grid[x][y+1][0] + gap_pen)
            if best == score + grid[x][y][0]:
                grid[x+1][y+1] = (best, 'd')
            elif best == grid[x+1][y][0] + gap_pen:
                grid[x+1][y+1] = (best, 'v')
            else:
                grid[x+1][y+1] = (best, 'h')

    for item in grid:
        out = ''
        for place in item:
            out += f'{place[0]:>8} '
        print(out)

    return grid


def global_backtrack(table, string1, string2):
    len1 = len(string1)
    len2 = len(string2)
    seq1 = ''
    seq2 = ''
    path = ''
    current = table[len2][len1]
    score = current[0]

    while current[1] != 's':
        if current[1] == 'd':
            seq1 = string1[len1-1] + seq1
            seq2 = string2[len2-1] + seq2
            len1 -= 1
            len2 -= 1
            path = 'D' + path
        elif current[1] == 'h':
            seq1 = '-' + seq1
            seq2 = string2[len2-1] + seq2
            len2 -= 1
            path = 'H' + path
        elif current[1] == 'v':
            seq1 = string1[len1 - 1] + seq1
            seq2 = '-' + seq2
            len1 -= 1
            path = 'V' + path
        current = table[len2][len1]

    print(f"Alignment score: {score}")
    print(f"Forward path is: {path}")

    print("Aligned sequences:")
    print(seq1)
    print(seq2)


def local_align(string1, string2, mat, o, gap_pen):
    grid = [[None for x in range(len(string1) + 1)] for y in range(len(string2) + 1)]
    grid[0][0] = (0, 's')

    for i in range(1, len(string1) + 1):
        grid[0][i] = (0, 'h')

    for i in range(1, len(string2) + 1):
        grid[i][0] = (0, 'v')

    for x in range(len(string2)):
        for y in range(len(string1)):
            score = get_score(mat, o, (string1[y], string2[x]))
            best = max(grid[x][y][0] + score, grid[x + 1][y][0] + gap_pen, grid[x][y + 1][0] + gap_pen)
            if best == score + grid[x][y][0]:
                if best < 0:
                    best = 0
                grid[x + 1][y + 1] = (best, 'd')
            elif best == grid[x + 1][y][0] + gap_pen:
                if best < 0:
                    best = 0
                grid[x + 1][y + 1] = (best, 'v')
            else:
                if best < 0:
                    best = 0
                grid[x + 1][y + 1] = (best, 'h')

    for item in grid:
        out = ''
        for place in item:
            out += f'{place[0]:>8} '
        print(out)

    return grid


def local_backtrack(table, string1, string2):
    max_n = 0
    for i in range(len(table)):
        for j in range(len(table[0])):
            if table[i][j][0] > max_n:
                max_n = table[i][j][0]
                placement = (i, j)

    len1 = placement[0]
    len2 = placement[1]
    seq1 = ''
    seq2 = ''
    path = ''
    current = table[placement[0]][placement[1]]
    score = current[0]

    while current[1] != 's':
        if current[0] == 0:
            break
        if current[1] == 'd':
            seq1 = string1[len1 - 1] + seq1
            seq2 = string2[len2 - 1] + seq2
            len1 -= 1
            len2 -= 1
            path = 'D' + path
        elif current[1] == 'h':
            seq1 = '-' + seq1
            seq2 = string2[len2 - 1] + seq2
            len2 -= 1
            path = 'H' + path
        elif current[1] == 'v':
            seq1 = string1[len1 - 1] + seq1
            seq2 = '-' + seq2
            len1 -= 1
            path = 'V' + path
        current = table[len2][len1]

    print(f"Alignment score: {score}")
    print(f"Forward path is: {path}")

    print("Aligned sequences:")
    print(seq1)
    print(seq2)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        # translation("ATATTGCAGTAGTAGCGAGTC")
        s1 = extract_string("Sample-Inputs/sequenceA1.txt")
        s2 = extract_string("Sample-Inputs/sequenceA2.txt")
        print("Sequences:")
        print(s1)
        print(s2)
        matrix, order = extract_matrix("Files/BLOSUM62.txt")
        gap = int(input("Enter the gap penalty: "))
        print("Starting...\n")
        grid = local_align(s1, s2, matrix, order, gap)
        local_backtrack(grid, s1, s2)
        # semiglobal_backtrack(grid, s1, s2)

        # grid = semiglobal_align(s1, s2, matrix, order, gap)
        # semiglobal_backtrack(grid, s1, s2)
        # grid = global_align(s1, s2, matrix, order, gap)
        # global_backtrack(grid, s1, s2)
        # print(get_score(matrix, order, ("A", "C")))
        # print(get_score(matrix, order, ("C", "A")))
        # print(get_score(matrix, order, ("Y", "E")))
        # print(s1)
        # print(s2)
    else:
        form = input("Will these inputs be nucleotide or peptide (n or p)? ").strip()[0]
        infile, infile2, mat_file = get_files()
        mut_prob = "Is the file a mutation probability PAM matrix (y or n)? "
        if mut_prob.lower() == "y":
            units = input("Please enter the units of divergence: ")
        s1 = extract_string(infile)
        s2 = extract_string(infile2)
        if form.lower() == "n":
            s1 = translation(s1)
            s2 = translation(s2)
            print("\n\n")
            print("Translated sequences:")
            print(s1)
            print(s2)
            print("\n\n")

        matrix, order = extract_matrix(mat_file)

        alignment = input("What type of alignment do you wish to use?\nGlobal, Local, Semi-Global (g, l, or s)? ").strip()[0]
        gap = int(input("Enter the gap penalty: "))
        print("Starting...\n")

        if alignment.lower() == "g":
            grid = global_align(s1, s2, matrix, order, gap)
            global_backtrack(grid, s1, s2)
        elif alignment.lower() == "s":
            grid = semiglobal_align(s1, s2, matrix, order, gap)
            semiglobal_backtrack(grid, s1, s2)
        elif alignment.lower() == "l":
            grid = local_align(s1, s2,matrix, order, gap)
            local_backtrack(grid, s1, s2)
        else:
            print("That is not a valid input.")
            exit()


