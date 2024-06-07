import os
import sys
import json
import pandas as pd
import logging
import time


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

    def createContigsInfo(self, virusKmerPool):
        # build contigsInfo to store: contig sequence, length of contig, virus kmers that exist in the contig (with location of kmer in contig)
        contigsInfo = []

        virusKmerPoolSet = set(virusKmerPool)
        for id, contig in enumerate(self.contigs):
            contigLen = len(contig)

            kmerCount = 0
            contigInfo = {"contig": contig, "length": contigLen, "v-kmers": []}
            contigSet = set(
                contig[i : i + self.k] for i in range(len(contig) - self.k + 1)
            )  # Create set of kmers in contig
            commonKmers = virusKmerPoolSet & contigSet  # Find common kmers
            if len(commonKmers) > 0:
                print(f"There are matched kmers in virus for contig: {id}")
                logging.info(f"\tThere are matched kmers in virus for contig: {id}")
            else:
                contigSpecificDistance = []
                localDistance = []
                for vkmer in virusKmerPool:
                    subLocalDistance = []
                    for ckmer in contigSet:
                        subLocalDistance.append(self.hammingDistance(vkmer, ckmer))
                    sldt = sum(subLocalDistance)
                    localDistance.append(sldt)
                contigSpecificDistance = sum(localDistance)
                print(f"Contig: {id} | Distance: {contigSpecificDistance}")

    def align(self, virusSequence):
        virusKmerPool = self.virusesToKmers(virusSequence)

        contigsInfo = self.createContigsInfo(virusKmerPool)

        possibleBestContigs = []
        mostVKmers = 0

        numGoodContigs = 0
        # find the contig with the most virus-kmers (count only)
        for contig in contigsInfo:
            if contig["kmerCount"] > mostVKmers:
                mostVKmers = contig["kmerCount"]
        print(f"mostVkmers: {mostVKmers}")

        # add all of the contigs that have mostVKmers count
        for contig in contigsInfo:
            if contig["kmerCount"] == mostVKmers:
                possibleBestContigs.append(contig)
            if contig["kmerCount"] == mostVKmers and mostVKmers >= 1:
                numGoodContigs += 1
        # print(f"possible contigs: {possibleBestContigs}")

        return numGoodContigs

    def search(self):
        logging.info("Search for Viruses: ")
        viruses = self.viruses
        for virus, virusData in viruses.items():
            print(f"Current Virus: {virus}")
            virusSequence = virusData["sequence"]
            start = time.time()
            numGoodContigs = self.align(virusSequence)
            end = time.time()
            logging.info(
                f"\tThere were {numGoodContigs} that had viral sequence data for virus: {virus}.\n\tIt took {end-start} seconds to digest this virus"
            )
