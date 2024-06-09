import json
import os
import logging
from multiprocessing import pool

# In your multiprocessing module
logger = logging.getLogger(__name__)
logger.setLevel(logging.CRITICAL)

handler = logging.FileHandler("multiprocessing.log")
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)

logger.addHandler(handler)


class SearchString:
    def __init__(self, viruses, readsKmerPool, contigs, k):
        self.viruses = viruses
        self.contigs = contigs
        self.readsKmerPool = readsKmerPool
        self.k = k
        self.maxHammingDistance = 2

    def virusToKmers(self, sequence):
        kmerPool = {}
        for index in range(len(sequence) - self.k + 1):
            kmer = sequence[index : index + self.k]
            if kmer not in kmerPool:
                kmerPool[kmer] = [index]
            else:
                kmerPool[kmer].append(index)
        return kmerPool

    def kmerPoolsToFile(self, virusKmerPool):
        with open("data/logs/r-kmerPool.json", "w") as file:
            json.dump(self.readsKmerPool, file)

        with open("data/logs/v-kmerPool.json", "w") as file:
            json.dump(virusKmerPool, file)

    def hammingDistance(self, virus, contig):
        distance = 0
        for index in range(len(virus)):
            if contig[index] != virus[index]:
                distance += 1
        return distance

    def createContigsInfo(self, virus, virusKmerPool):
        contigsInfo = []
        for id, contig in enumerate(self.contigs):
            contigLen = len(contig)
            contigKmers = {
                contig[i : i + self.k] for i in range(contigLen - self.k + 1)
            }
            kmerCount = 0
            contigInfo = {
                "virus": virus,
                "contigId": id + 1,
                "contig": contig,
                "length": contigLen,
                "v-kmers": [],
            }
            for virusKmer in virusKmerPool:
                for contigKmer in contigKmers:
                    distance = self.hammingDistance(virusKmer, contigKmer)
                    if (
                        distance <= self.maxHammingDistance
                    ):  # max_distance is the maximum allowed Hamming distance
                        kmerCount += 1
                        contigInfo["v-kmers"].append(
                            {
                                virusKmer: {
                                    "indexOfVKmerInContig": [
                                        contig.index(contigKmer),
                                        contig.index(contigKmer) + self.k,
                                    ],
                                    "hammingDistance": distance,
                                }
                            }
                        )
            """for kmer in virusKmerPool:
                if kmer in contigKmers:
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
                    )"""
            contigInfo["kmerCount"] = kmerCount
            contigsInfo.append(contigInfo)
        return contigsInfo

    def searchString(self):
        logging.info("Search For Viruses:\n")
        virusesInBiosample = []
        for virus_id, virus in self.viruses.items():
            contigsExistInVirus = []
            virusKmerPool = self.virusToKmers(virus["sequence"])
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
                f"\t{len(contigsExistInVirus)} contigs aligns with {virus['name']} "
            )
        with open("data/output_data/virusesInBiosample.json", "w") as file:
            json.dump(virusesInBiosample, file)
        return virusesInBiosample
