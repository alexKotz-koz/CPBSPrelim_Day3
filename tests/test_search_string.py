import os
import sys
import json
import unittest
import pandas as pd

# os management for testing
sys.path.insert(0, "../src/components")

from searchForVirsues_old import SearchString


class TestSearchString(unittest.TestCase):
    def setUp(self):
        # test data for testQueryToKmers
        rawQueryData = {
            "id": ["INITIAL_QUERY"],
            "sequence": ["GGATC"],
        }
        rawQueryData["length"] = [len(seq) for seq in rawQueryData["sequence"]]
        self.queryDataForKmers = pd.DataFrame(rawQueryData)

        # test data for testCreateContigInfo and testAlign
        queryData = []
        with open("../src/data/QUERY.fasta", "r") as inputQueryFile:
            query = inputQueryFile.readlines()
            for index, line in enumerate(query):
                if ">" == line[0]:
                    queryString = query[index + 1]
                    queryData.append(
                        {
                            "id": line.lstrip(">").rstrip("\n"),
                            "sequence": queryString.rstrip("\n"),
                        }
                    )
        self.dfQueryData = pd.DataFrame(queryData)
        self.dfQueryData["length"] = self.dfQueryData["sequence"].str.len()

        self.contigs = [
            "ACTGGATCTTCAG",
            "ACTGGATCTTCG",
            "TTTAGCCTTCG",
            "TTTAGCCTTCAG",
            "TTTAGCCTTATC",
        ]

        self.readsKmerPool = {
            "ACTGG": {"Read1": [{0: 5}]},
            "CTGGA": {"Read1": [{1: 6}]},
            "TGGAT": {"Read1": [{2: 7}]},
            "GGATC": {"Read1": [{3: 8}]},
            "GATCT": {"Read1": [{4: 9}]},
            "ATCTT": {"Read1": [{5: 10}]},
            "TCTTC": {"Read1": [{6: 11}]},
            "CTTCA": {"Read1": [{7: 12}]},
            "TTCAG": {"Read1": [{8: 13}]},
            "CTAGC": {"Read2": [{0: 5}]},
            "TAGCC": {"Read2": [{1: 6}]},
            "AGCCT": {"Read2": [{2: 7}], "Read3": [{0: 5}]},
            "GCCTT": {"Read2": [{3: 8}], "Read3": [{1: 6}]},
            "CCTTA": {"Read2": [{4: 9}]},
            "CTTAT": {"Read2": [{5: 10}]},
            "TTATC": {"Read2": [{6: 11}]},
            "CCTTC": {"Read3": [{2: 7}]},
            "CTTCG": {"Read3": [{3: 8}]},
            "TTTAG": {"Read4": [{0: 5}]},
            "TTAGC": {"Read4": [{1: 6}]},
            "TAGCT": {"Read4": [{2: 7}]},
            "AGCTA": {"Read4": [{3: 8}]},
            "GCTAG": {"Read4": [{4: 9}]},
        }

        self.contigsInfoTest = [
            {
                "contig": "ACTGGATCTTCAG",
                "length": 13,
                "q-kmers": [
                    {"GGATC": {"index": [3, 8]}},
                    {"TGGAT": {"index": [2, 7]}},
                    {"ACTGG": {"index": [0, 5]}},
                    {"GATCT": {"index": [4, 9]}},
                    {"CTGGA": {"index": [1, 6]}},
                ],
                "kmerCount": 5,
            },
            {
                "contig": "ACTGGATCTTCG",
                "length": 12,
                "q-kmers": [
                    {"GGATC": {"index": [3, 8]}},
                    {"TGGAT": {"index": [2, 7]}},
                    {"ACTGG": {"index": [0, 5]}},
                    {"GATCT": {"index": [4, 9]}},
                    {"CTGGA": {"index": [1, 6]}},
                ],
                "kmerCount": 5,
            },
            {"contig": "TTTAGCCTTCG", "length": 11, "q-kmers": [], "kmerCount": 0},
            {"contig": "TTTAGCCTTCAG", "length": 12, "q-kmers": [], "kmerCount": 0},
            {"contig": "TTTAGCCTTATC", "length": 12, "q-kmers": [], "kmerCount": 0},
        ]

        self.readsInContig = [
            {"ACTGG": {"Read1": [{0: 5}]}},
            {"CTGGA": {"Read1": [{1: 6}]}},
            {"TGGAT": {"Read1": [{2: 7}]}},
            {"GGATC": {"Read1": [{3: 8}]}},
            {"GATCT": {"Read1": [{4: 9}]}},
            {"ATCTT": {"Read1": [{5: 10}]}},
            {"TCTTC": {"Read1": [{6: 11}]}},
            {"CTTCA": {"Read1": [{7: 12}]}},
            {"TTCAG": {"Read1": [{8: 13}]}},
        ]

        self.k = 5

        self.searchStringInstance = SearchString(
            self.dfQueryData, self.readsKmerPool, self.contigs, self.k
        )

    def testQueryToKmers(self):
        # secondary searchString instance for smaller query to kmer test
        self.searchStringInstance2 = SearchString(
            self.queryDataForKmers, self.readsKmerPool, self.contigs, self.k
        )
        kmerPool = self.searchStringInstance2.queryToKmers()
        kmerPoolTest = {"GGATC": [{0: 5}]}
        self.assertDictEqual(kmerPoolTest, kmerPool)

    def testCreateContigInfo(self):
        with open("test_data/q-kmerPool_test.json", "r") as file:
            qKmerPool = json.load(file)

        contigsInfo = self.searchStringInstance.createContigsInfo(qKmerPool)
        self.assertEqual(contigsInfo, self.contigsInfoTest)


if __name__ == "__main__":
    unittest.main()
