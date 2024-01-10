from Bio.Data import IUPACData
from itertools import product, permutations


def replace_slash(seq: str) -> str:
    """Replace slash notation with IUPAC ambiguous DNA code."""
    def invert(d):
        result = {}
        for key, value in d.copy().items():
            if len(value) == 1:
                result[value] = key
                continue
            for permutation in permutations(value):
                result[f'({"/".join(permutation)})'] = key
        return result

    for key, value in invert(IUPACData.ambiguous_dna_values).items():
        seq = seq.replace(key, value)
    return seq


def disambiguate_sequence(seq: str) -> list[str]:
    """Given an ambiguous DNA sequence, return a list of all possible sequences."""
    seq = replace_slash(seq)
    return list(map("".join, product(*map(IUPACData.ambiguous_dna_values.get, seq))))
