import sys
import os

# required os set up for testing
currentDir = os.path.dirname(os.path.abspath(__file__))
parentDir = os.path.dirname(currentDir)
sys.path.insert(0, parentDir)

import json
import pandas as pd


class SearchString:
    def __init__(self, virus, readsKmerPool, contigs, k):
        self.virus = virus
        self.contigs = contigs
        self.readsKmerPool = readsKmerPool
        self.k = 10
        # error checking for test cases, where toy dataset reads and kmers are smaller than the original reads data
        if k < self.k:
            self.k = k

    def virusToKmers(self):
        kmerPool = {}
        sequence = self.virus["sequence"]
        for index, base in enumerate(sequence):
            kmer = sequence[index : index + self.k]

            if len(kmer) >= self.k:
                if kmer not in kmerPool:
                    kmerPool[kmer] = [{index: index + self.k}]
                else:
                    kmerPool[kmer].append({index: index + self.k})
        return kmerPool

    def kmerPoolsToFile(self, virusKmerPool):

        # write reads and query kmerpool to files for analysis
        with open("data/logs/v-kmerPool.json", "w") as file:
            json.dump(virusKmerPool, file)

        with open("data/logs/r-kmerPool.json", "w") as file:
            json.dump(self.readsKmerPool, file)

    def createContigsInfo(self, virusKmerPool):
        # build contigsInfo to store: contig sequence, length of contig, virus string kmers that exist in the contig (with location of kmer in contig)
        contigsInfo = []
        for contig in self.contigs:
            contigLen = len(contig)
            kmerCount = 0
            contigInfo = {"contig": contig, "length": contigLen, "v-kmers": []}
            for index, kmer in enumerate(virusKmerPool):
                if kmer in contig and kmer not in contigInfo["v-kmers"]:
                    kmerCount += 1
                    contigInfo["v-kmers"].append(
                        {
                            kmer: {
                                "index": [
                                    contig.index(kmer),
                                    contig.index(kmer) + self.k,
                                ]
                            }
                        }
                    )
            contigInfo["kmerCount"] = kmerCount
            contigsInfo.append(contigInfo)
        return contigsInfo

    def align(self):
        # break query string into kmers
        virusKmerPool = self.virusToKmers()
        # build contigsInfo (see method for description)
        contigsInfo = self.createContigsInfo(virusKmerPool)
        # find longest contig that contains the most query string kmers
        possibleBestContigs = []
        mostVKmers = 0
        longestContigMostVKmers = ""
        # find mostVKmers
        for contig in contigsInfo:
            if contig["kmerCount"] > mostVKmers:
                mostVKmers = contig["kmerCount"]

        for contig in contigsInfo:
            if contig["kmerCount"] == mostVKmers:
                possibleBestContigs.append(contig)

        chronologicalOrder = {}
        for contig in possibleBestContigs:
            for vKmer in contig["v-kmers"]:
                startIndex = list(vKmer.values())[0]["index"][0]
                if contig["contig"] not in chronologicalOrder:
                    chronologicalOrder[contig["contig"]] = [startIndex]
                else:
                    chronologicalOrder[contig["contig"]].append(startIndex)

        if len(chronologicalOrder) == 1:  # passed
            longestContigMostVKmers = list(chronologicalOrder.keys())[0]
        else:
            # minimum difference between starting indecies
            minDifference = float("inf")
            minContigs = []
            maxKmerCount = float("-inf")
            maxContig = None
            contigsInfoDict = {item["contig"]: item for item in contigsInfo}
            for contig, indices in chronologicalOrder.items():
                # Calculate the differences between consecutive elements
                diffs = [j - i for i, j in zip(indices[:-1], indices[1:])]
                # if the max difference is smaller than the current minimum, update tracking variables
                if max(diffs) < minDifference:
                    minDifference = max(diffs)
                    minContigs = [contig]

                # if the max difference is equal to the current minimum, add the contig to minContigs
                elif max(diffs) == minDifference:
                    minContigs.append(contig)

            # if there is a tie, find the longest contig in contigsInfo
            if len(minContigs) > 1:
                maxKmerCount = float("-inf")
                for contig in minContigs:
                    if contigsInfoDict[contig]["length"] > maxKmerCount:
                        maxKmerCount = contigsInfoDict[contig]["length"]
                        maxContig = contig
            longestContigMostVKmers = maxContig
            # index of contig in self.contigs returned from contigCreation
            contigIndex = self.contigs.index(longestContigMostVKmers)

        rootDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        filePathALLELES = os.path.join(rootDir, "data/output/ALLELES.fasta")
        with open(filePathALLELES, "w") as file:
            file.write(longestContigMostVKmers)

        # find the reads that exist in each contig for ALLELES.fasta and output.aln
        readsInContig = []
        for read in self.readsKmerPool:
            if read in longestContigMostVKmers:
                readsInContig.append({read: self.readsKmerPool[read]})

        filePathReadsInContig = os.path.join(rootDir, "data/logs/readsInContig.json")
        with open(filePathReadsInContig, "w") as file:
            json.dump(readsInContig, file)

        # find the query k-mers that exist in the final contig and where they exist
        vKmerInFinalContig = []
        for vKmer in virusKmerPool:
            if vKmer in longestContigMostVKmers:
                vKmerInFinalContig.append({vKmer: virusKmerPool[vKmer]})
        filePathVKmersInContig = os.path.join(rootDir, "data/logs/vKmerInContig.json")
        with open(filePathVKmersInContig, "w") as file:
            json.dump(vKmerInFinalContig, file)

        return (
            contigsInfo,
            longestContigMostVKmers,
            contigIndex,
            readsInContig,
            vKmerInFinalContig,
        )

    def searchString(self):
        virusKmerPool = self.virusToKmers()
        self.kmerPoolsToFile(virusKmerPool=virusKmerPool)
        (
            contigsInfo,
            longestContigMostVKmers,
            contigIndex,
            readsInContig,
            vKmerInFinalContig,
        ) = self.align()
        return (
            contigsInfo,
            longestContigMostVKmers,
            contigIndex,
            readsInContig,
            vKmerInFinalContig,
        )
