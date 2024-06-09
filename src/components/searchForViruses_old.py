import json
import os


class SearchString:
    def __init__(self, viruses, readsKmerPool, contigs, k):
        self.viruses = viruses
        self.contigs = contigs
        self.readsKmerPool = readsKmerPool
        self.k = k

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

    def createContigsInfo(self, virus, virusKmerPool):
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
            for kmer in virusKmerPool:
                if kmer in contig:
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
        with open("data/output_data/virusesInBiosample.json", "w") as file:
            json.dump(virusesInBiosample, file)
        return virusesInBiosample
