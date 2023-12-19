from typing import List

from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools


class Target:

    def __init__(
        self,
        sequence: Seq,
        pam: Seq,
        position: int,
        strand: str = "fwd",
        target_length: int = 20
    ):
        self.sequence = sequence
        self.pam = pam
        self.position = position  # 0-based index
        self.strand = strand
        self.target_length = target_length

    def __len__(self):
        return len(self.pam) + len(self.sequence)

    def __str__(self):
        data = self.json()
        return f"{data['position']} {data['targetSequence']} {data['targetPam']}"

    def __lt__(self, other):
        """Sort by strand and then by position"""  # We're probably going to implement more sorting options later, so this function is basically a placeholder.
        if self.strand == other.strand:
            return self.position < other.position
        return self.strand == "fwd"  # reverse strand comes after forward

    def json(self):
        return {
            "position": self.position + 1,  # for 1-based positioning
            "strand": self.strand,
            "sequence": str(self.sequence) if self.strand == "fwd" else str(self.sequence.reverse_complement()),
            "targetSequence": str(self.sequence),
            "pam": str(self.pam) if self.strand == "fwd" else str(self.pam.reverse_complement()),
            "targetPam": str(self.pam),  # PAM as it is in the target sequence
            "targetLength": self.target_length
        }


class TargetList:

        def __init__(self, sequence_record: SeqRecord, target_length: int = 20, pam: str = "NGG"):
            self.name = sequence_record.name
            self.description = sequence_record.description
            self.sequence = sequence_record.seq
            self.pam = pam
            self.target_length = target_length
            self.targets = []
            self._find_targets()

        def append(self, target: Target):
            if len(target) == self.target_length + len(self.pam):
                self.targets.append(target)

        def __len__(self):
            return len(self.targets)

        def _replace_wildcards(self, pam: str) -> List[Seq]:
            # If there are no more wildcards, return the sequence
            if 'N' not in pam:
                return [Seq(pam)]

            # Replace the first wildcard with each possible nucleotide and recursively call the function
            possible_sequences = []
            for nucleotide in ['A', 'T', 'C', 'G']:
                new_sequences = self._replace_wildcards(pam.replace('N', nucleotide, 1))
                possible_sequences.extend(new_sequences)

            return possible_sequences

        def _cut_sequence(self, motif_position: int, strand: str = "fwd") -> Seq:
            start_position = max(0, motif_position - self.target_length)
            end_position = motif_position

            if strand == "rev":
                start_position = max(0, (motif_position + len(self.pam)))
                end_position = motif_position + self.target_length + len(self.pam)

            return self.sequence[start_position:end_position]

        def _find_targets(self):
            for strand, pam in itertools.product(["fwd", "rev"], self._replace_wildcards(self.pam)):
                pam = Seq(pam) if strand == "fwd" else Seq(pam).reverse_complement()
                motif = motifs.create([pam])
                for motif_instance in motif.instances.search(self.sequence):
                    self.append(
                        Target(
                            sequence=self._cut_sequence(motif_instance[0], strand),
                            pam=pam,
                            target_length=self.target_length,
                            position=motif_instance[0],
                            strand=strand
                        )
                    )

        def json(self):
            return {
                "name": self.name,
                "description": self.description,
                "sequence": str(self.sequence),
                "pam": self.pam,
                "targetLength": self.target_length,
                "targets": [target.json() for target in sorted(self.targets)]
            }
