import sys
import os
import logging

# required os set up for testing
currentDir = os.path.dirname(os.path.abspath(__file__))
parentDir = os.path.dirname(currentDir)
sys.path.insert(0, parentDir)

import json
import pandas as pd


class SearchString:
    def __init__(self, viruses, readsKmerPool, contigs, k):
        self.viruses = viruses
        self.contigs = contigs
        self.readsKmerPool = readsKmerPool
        self.k = k

    def virusToKmers(self, sequence):
        kmerPool = {}
        for index, base in enumerate(sequence):
            kmer = sequence[index : index + self.k]
            if len(kmer) >= self.k:
                if kmer not in kmerPool:
                    kmerPool[kmer] = [{index: index + self.k}]
                else:
                    kmerPool[kmer].append({index: index + self.k})
        return kmerPool

    def kmerPoolsToFile(self, virusKmerPool):

        with open("data/logs/r-kmerPool.json", "w") as file:
            json.dump(self.readsKmerPool, file)

        # write reads and query kmerpool to files for analysis
        with open("data/logs/v-kmerPool.json", "w") as file:
            json.dump(virusKmerPool, file)

    def createContigsInfo(self, virus, virusKmerPool):
        # build contigsInfo to store: contig sequence, length of contig, query string kmers that exist in the contig (with location of kmer in contig)
        contigsInfo = []
        for id, contig in enumerate(self.contigs):
            contigLen = len(contig)
            kmerCount = 0
            contigInfo = {
                "virus": virus,
                "contigId": id + 1,
                "contig": contig,
                "length": contigLen,
                "v-kmers": [],
            }
            for index, kmer in enumerate(virusKmerPool):
                if kmer in contig and kmer not in contigInfo["v-kmers"]:
                    kmerCount += 1
                    contigInfo["v-kmers"].append(
                        {
                            kmer: {
                                "indexOfVKmerInContig": [
                                    contig.index(kmer),
                                    contig.index(kmer) + self.k,
                                ]
                            }
                        }
                    )
            contigInfo["kmerCount"] = kmerCount
            contigsInfo.append(contigInfo)
        return contigsInfo

    def searchString(self):
        virusesInBiosample = []
        for ssr, virus in self.viruses.items():
            contigsExistInVirus = []
            virusKmerPool = self.virusToKmers(virus["sequence"])
            # build contigsInfo (see method for description)
            contigsInfo = self.createContigsInfo(virus["name"], virusKmerPool)
            for contig in contigsInfo:
                if contig["kmerCount"] > 0:
                    contigsExistInVirus.append(contig)
            virusesInBiosample.append(
                {
                    "virus": virus["name"],
                    "numContigsInVirus": len(contigsExistInVirus),
                    "contigsInVirus": contigsExistInVirus,
                }
            )
            logging.info(
                f"There are {len(contigsExistInVirus)} contigs that align with virus: {virus['name']} of length {len(virus['sequence'])}"
            )
        with open("data/output_data/virusesInBiosample.json", "w") as file:
            json.dump(virusesInBiosample, file)
        return virusesInBiosample
