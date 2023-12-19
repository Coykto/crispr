from typing import List

from Bio.Seq import Seq
from Bio import motifs
from Bio.SeqIO import parse


class Target:

    def __init__(
        self,
        sequence: Seq,
        pam: Seq,
        target_length: int,
        position: int,
        strand: str = 'fwd'
    ):
        self.strand = strand
        self.sequence = sequence
        self.pam = pam
        self.target_length = target_length
        self.position = position + 1

    def __len__(self):
        return len(self.pam) + len(self.sequence)

    def __str__(self):
        return f"{self.position} {self.sequence} {self.pam}"


class TargetList():

    def __init__(self, target_length: int, pam: str):
        self.targets = []
        self.target_length = target_length
        self.pam = pam

    def append(self, target: Target):
        if len(target) != self.target_length + len(self.pam):
            return
        self.targets.append(target)

    def __len__(self):
        return len(self.targets)

    def __str__(self):
        return '\n'.join([str(target) for target in sorted(self.targets, key=lambda target: target.position)])


def find_targets(sequence: Seq, pam: str, target_length: int) -> TargetList:

    possible_pams = replace_wildcards(pam)

    targets = TargetList(target_length, pam)
    for pam in possible_pams:
        motif = motifs.create([pam])
        for motif_instance in motif.instances.search(sequence):
            motif_position = motif_instance[0]
            targets.append(
                Target(
                    sequence=sequence[max(0, motif_position - target_length):motif_position],
                    pam=pam,
                    target_length=target_length,
                    position=motif_position,
                    strand="fwd"
                )
            )
    return targets


def find_reverse_targets(sequence: Seq, pam: str, target_length: int) -> TargetList:
    reverse_pam = str(Seq(pam).reverse_complement())

    # possible_pams = replace_wildcards(reverse_pam)
    possible_pams = [Seq("CCT")]

    targets = TargetList(target_length, pam)
    for pam in possible_pams:
        motif = motifs.create([pam])
        for motif_instance in motif.instances.search(sequence):
            motif_position = motif_instance[0]

            primer_start_position = max(0, (motif_position + len(pam)))
            primer_end_position = motif_position + target_length + len(pam)

            targets.append(
                Target(
                    sequence=sequence[primer_start_position:primer_end_position],
                    pam=pam,
                    target_length=target_length,
                    position=motif_position,
                    strand="rev"
                )
            )
    return targets


def replace_wildcards(pam: str) -> List[Seq]:
    # If there are no more wildcards, return the sequence
    if 'N' not in pam:
        return [Seq(pam)]

    # Replace the first wildcard with each possible nucleotide and recursively call the function
    possible_sequences = []
    for nucleotide in ['A', 'T', 'C', 'G']:
        new_sequences = replace_wildcards(pam.replace('N', nucleotide, 1))
        possible_sequences.extend(new_sequences)

    return possible_sequences



if __name__ == '__main__':
    PAM = "NGG"
    sequence = Seq("CTTCCTTTGTCCCCAATCTGGGCGCGCGCCGGCGCCCCCTGGCGGCCTAAGGACTCGGCGCGCCGGAAGTGGCCAGGGCGGGGGCGACCTCGGCTCACAGCGCGCCCGGCTATTCTCGCAGCTCACCATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCC")
    targets = find_targets(sequence, PAM, 20)
    print("FWD:")
    print(targets)

    print("REV:")
    print(find_reverse_targets(sequence, PAM, 20))


    '''
    fwd
    1 TCCCCAATCTGGGCGCGCGC CGG
    2 GGGCGCGCGCCGGCGCCCCC TGG
    3 CGCGCGCCGGCGCCCCCTGG CGG
    4 CGGCGCCCCCTGGCGGCCTA AGG
    5 CCCCTGGCGGCCTAAGGACT CGG
    6 GGCCTAAGGACTCGGCGCGC CGG
    7 AGGACTCGGCGCGCCGGAAG TGG
    8 TCGGCGCGCCGGAAGTGGCC AGG
    9 CGGCGCGCCGGAAGTGGCCA GGG
    10 CGCGCCGGAAGTGGCCAGGG CGG
    11 GCGCCGGAAGTGGCCAGGGC GGG
    12 CGCCGGAAGTGGCCAGGGCG GGG
    13 GCCGGAAGTGGCCAGGGCGG GGG
    14 GGCCAGGGCGGGGGCGACCT CGG
    15 ACCTCGGCTCACAGCGCGCC CGG
    16 GCTATTCTCGCAGCTCACCA TGG
    17 GCCGCGCTCGTCGTCGACAA CGG
    18 CTCGTCGTCGACAACGGCTC CGG
    19 CAACGGCTCCGGCATGTGCA AGG
    20 GGCTCCGGCATGTGCAAGGC CGG
    21 CATGTGCAAGGCCGGCTTCG CGG
    22 ATGTGCAAGGCCGGCTTCGC GGG
    23 TCGCGGGCGACGATGCCCCC CGG
    24 CGCGGGCGACGATGCCCCCC GGG
    25 GGCCGTCTTCCCCTCCATCG TGG
    26 GCCGTCTTCCCCTCCATCGT GGG
    27 CCGTCTTCCCCTCCATCGTG GGG
    rev
    1 CGCGCCCAGATTGGGGACAA AGG
    2 GCCGGCGCGCGCCCAGATTG GGG
    3 CGCCGGCGCGCGCCCAGATT GGG
    4 GCGCCGGCGCGCGCCCAGAT TGG
    5 CTTAGGCCGCCAGGGGGCGC CGG
    6 CGAGTCCTTAGGCCGCCAGG GGG
    7 CCGAGTCCTTAGGCCGCCAG GGG
    8 GCCGAGTCCTTAGGCCGCCA GGG
    9 CGCCGAGTCCTTAGGCCGCC AGG
    10 TTCCGGCGCGCCGAGTCCTT AGG
    11 GCCCCCGCCCTGGCCACTTC CGG
    12 AGCCGAGGTCGCCCCCGCCC TGG
    13 GCCGGGCGCGCTGTGAGCCG AGG
    14 GGTGAGCTGCGAGAATAGCC GGG
    15 TGGTGAGCTGCGAGAATAGC CGG
    16 CGCGGCGATATCATCATCCA TGG
    17 GCCGTTGTCGACGACGAGCG CGG
    18 GAAGCCGGCCTTGCACATGC CGG
    19 GGCATCGTCGCCCGCGAAGC CGG
    20 GGAGGGGAAGACGGCCCGGG GGG
    21 TGGAGGGGAAGACGGCCCGG GGG
    22 ATGGAGGGGAAGACGGCCCG GGG
    23 GATGGAGGGGAAGACGGCCC GGG
    24 CGATGGAGGGGAAGACGGCC CGG
    25 CCCCACGATGGAGGGGAAGA CGG
    
    
    
    Seq('GCCCAGATTGGGGACAAAGG')
    '''