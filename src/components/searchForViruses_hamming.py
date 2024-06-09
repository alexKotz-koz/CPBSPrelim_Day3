import os
import sys
import json
import pandas as pd
import logging
import time
from utils.utils import ToKmers


class SearchForViruses:
    def __init__(self, viruses, readsKmerPool, contigs, k):
        self.viruses = viruses
        self.readsKmerPool = readsKmerPool
        self.contigs = contigs
        self.k = k

    def virusesToKmers(self, virusSequence):
        kmerPool = {}
        for key, virus in self.viruses.items():
            sequence = virus["sequence"]

        for index, base in enumerate(sequence):
            kmer = sequence[index : index + self.k]
            if len(kmer) >= self.k:
                if kmer not in kmerPool:
                    kmerPool[kmer] = [{index: index + self.k}]
                else:
                    kmerPool[kmer].append({index: index + self.k})
        return kmerPool

    def hammingDistance(self, virus, contig):
        distance = 0
        for index in range(len(virus)):
            if contig[index] != virus[index]:
                distance += 1
        return distance

    def align(self, virusKmerPool):
        virusKmerPoolSet = set(virusKmerPool)
        for id, contig in enumerate(self.contigs):
            contigLen = len(contig)
            kmerCount = 0
            contigInfo = {"contig": contig, "length": contigLen, "v-kmers": []}
            contigSet = set(
                contig[i : i + self.k] for i in range(len(contig) - self.k + 1)
            )  # Create set of kmers in contig

    def search(self):
        logging.info("Search for Viruses: ")
        viruses = self.viruses
        for virus, virusData in viruses.items():
            virusSequence = virusData["sequence"]
            start = time.time()
            numGoodContigs = self.align(virusSequence)
            end = time.time()
