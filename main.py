import re

from itertools import product, permutations


def inverted(d):
    result = {}
    for key, value in d.copy().items():
        if len(value) == 1:
            result[value] = key
            continue
        for permutation in permutations(value):
            result[f'({"/".join(permutation)})'] = key
    return result

def replace_slash(seq):
    inverted_ambiguous_dna_values = inverted(IUPACData.ambiguous_dna_values)
    for key, value in inverted_ambiguous_dna_values.items():
        seq = seq.replace(key, value)
    return seq



if __name__ == '__main__':
    from Bio.Data import IUPACData

    seq = "N(A/G)NG(A/G)T(G/A)"
    seq = replace_slash(seq)
    # seq = seq.replace("(A/G)", "R")
    l = list(map("".join, product(*map(IUPACData.ambiguous_dna_values.get, seq))))
    # awd = 23


    a = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "M": "AC"
    }
    b = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "(A/C)": "M",
        "(C/A)": "M",
    }
    g = inverted(IUPACData.ambiguous_dna_values)
    awd =23

