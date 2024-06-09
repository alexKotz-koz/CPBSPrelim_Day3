import os
import sys
import json
import pandas as pd
import logging
import time


class SearchForVirusesHamming:
    def __init__(self, viruses, reads, k):
        self.viruses = viruses
        self.reads = reads
        self.k = k

    def virusesToKmers(self, k, sequence):
        kmerPool = {}
        print(k)
        for index, base in enumerate(sequence):
            kmer = sequence[index : index + k]
            if len(kmer) >= k:
                if kmer not in kmerPool:
                    kmerPool[kmer] = [{index: index + k}]
                else:
                    kmerPool[kmer].append({index: index + k})
        return kmerPool

    def hammingDistance(self, virus, contig):
        distance = 0
        for index in range(len(virus)):
            if contig[index] != virus[index]:
                distance += 1
        return distance

    def align(self, virusKmerPool, read):
        contigKmers = {}
        for kmer in virusKmerPool:
            print(len(kmer))
            print(len(read))

    def search(self):
        logging.info("Search for Viruses: ")
        viruses = self.viruses
        reads = self.reads

        for virus, virusData in viruses.items():
            virusSequence = virusData["sequence"]
            print(virus)
            for index, read in reads.iterrows():
                readLen = len(read["sequence"])
                if len(virusSequence) < readLen:
                    print(f"Skipping virus: {virus}")
                    break
                else:
                    vKmers = self.virusesToKmers(k=readLen - 1, sequence=virusSequence)
                    self.align(vKmers, read["sequence"])
                    break
