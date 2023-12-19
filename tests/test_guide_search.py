import unittest
from Bio.Seq import Seq
from guides_search import Target, TargetList
from Bio.SeqRecord import SeqRecord

class TestTarget(unittest.TestCase):
    def setUp(self):
        self.target = Target(sequence=Seq("ATCG"), pam=Seq("NGG"), position=0)

    def test_init(self):
        self.assertEqual(str(self.target.sequence), "ATCG")
        self.assertEqual(str(self.target.pam), "NGG")
        self.assertEqual(self.target.position, 0)

    def test_len(self):
        self.assertEqual(len(self.target), 6)

    def test_str(self):
        self.assertEqual(str(self.target), "0 ATCG NGG")

    def test_lt(self):
        other = Target(sequence=Seq("ATCG"), pam=Seq("NGG"), position=1)
        self.assertTrue(self.target < other)

    def test_json(self):
        self.assertEqual(self.target.json(), {
            "position": 1,
            "strand": "fwd",
            "sequence": "ATCG",
            "targetSequence": "ATCG",
            "pam": "NGG",
            "targetPam": "NGG",
            "targetLength": 20
        })

class TestTargetList(unittest.TestCase):
    def setUp(self):
        self.target_list = TargetList(sequence_record=SeqRecord(Seq("ATCG")), target_length=4, pam="NGG")

    def test_init(self):
        self.assertEqual(str(self.target_list.sequence), "ATCG")
        self.assertEqual(self.target_list.pam, "NGG")
        self.assertEqual(self.target_list.target_length, 4)

    def test_append(self):
        target = Target(sequence=Seq("ATCG"), pam=Seq("NGG"), position=0)
        self.target_list.append(target)
        self.assertEqual(len(self.target_list), 1)

    def test_len(self):
        self.assertEqual(len(self.target_list), 0)

    def test_replace_wildcards(self):
        self.assertEqual(self.target_list._replace_wildcards("N"), [Seq("A"), Seq("T"), Seq("C"), Seq("G")])
        self.assertEqual(
            self.target_list._replace_wildcards("NAN"),
            [
                Seq('AAA'), Seq('AAT'), Seq('AAC'), Seq('AAG'), Seq('TAA'), Seq('TAT'), Seq('TAC'), Seq('TAG'),
                Seq('CAA'), Seq('CAT'), Seq('CAC'), Seq('CAG'), Seq('GAA'), Seq('GAT'), Seq('GAC'), Seq('GAG')
            ])

    def test_cut_sequence(self):
        self.assertEqual(self.target_list._cut_sequence(0), Seq(""))
        self.assertEqual(self.target_list._cut_sequence(4), Seq("ATCG"))

    def test_find_targets(self):
        self.target_list._find_targets()
        self.assertEqual(len(self.target_list), 2)

    def test_json(self):
        self.assertEqual(self.target_list.json(), {
            "name": "<unknown name>",
            "description": "<unknown description>",
            "sequence": "ATCG",
            "pam": "NGG",
            "targetLength": 4,
            "targets": []
        })

if __name__ == '__main__':
    unittest.main()