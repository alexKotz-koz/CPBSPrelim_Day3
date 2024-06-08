import os
import sys
import json
import pandas as pd
import numpy as np
import logging
import time

from concurrent.futures import ProcessPoolExecutor


class SearchForViruses:
    def __init__(self, viruses, contigs, k):
        self.viruses = viruses
        self.contigs = contigs
        self.k = k
        self.logDataDir = "./data/logs"
        self.outputDataDir = "./data/output_data"

    def smith_waterman(
        self, contig, virus, match_score=3, mismatch_score=-1, gap_penalty=-2
    ):
        """
        Perform Smith-Waterman alignment on two sequences.

        Args:
        - contig (str): First sequence to align.
        - virus (str): Second sequence to align.
        - match_score (int): Score for a match.
        - mismatch_score (int): Score for a mismatch.
        - gap_penalty (int): Penalty for opening a gap.

        Returns:
        - alignment_score (int): Score of the best alignment.
        - aligned_contig (str): First sequence with gaps for alignment.
        - aligned_virus (str): Second sequence with gaps for alignment.
        """
        # Initialize scoring matrix
        rows = len(contig) + 1
        cols = len(virus) + 1
        score_matrix = np.zeros((rows, cols), dtype=int)

        # Initialize traceback matrix
        traceback_matrix = np.zeros((rows, cols), dtype=int)

        # Starting from the second row and second column,
        for i in range(1, rows):
            for j in range(1, cols):
                match = score_matrix[i - 1][j - 1] + (
                    match_score if contig[i - 1] == virus[j - 1] else mismatch_score
                )
                delete = score_matrix[i - 1][j] + gap_penalty
                insert = score_matrix[i][j - 1] + gap_penalty
                score_matrix[i][j] = max(match, delete, insert, 0)
                traceback_matrix[i][j] = np.argmax([0, match, delete, insert])

        # Find the maximum score in the matrix
        max_score = np.max(score_matrix)

        # Find the index of the maximum score
        max_index = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)

        # Traceback to reconstruct the aligned sequences
        aligned_contig = ""
        aligned_virus = ""
        i, j = max_index
        startPosition = j
        while i > 0 and j > 0 and score_matrix[i][j] != 0:
            if traceback_matrix[i][j] == 1:
                aligned_contig = contig[i - 1] + aligned_contig
                aligned_virus = virus[j - 1] + aligned_virus
                i -= 1
                j -= 1
            elif traceback_matrix[i][j] == 2:
                aligned_contig = contig[i - 1] + aligned_contig
                aligned_virus = "-" + aligned_virus
                i -= 1
            else:
                aligned_contig = "-" + aligned_contig
                aligned_virus = virus[j - 1] + aligned_virus
                j -= 1
        endPosition = j
        # print(f"Done: {max_score}")
        return max_score, aligned_contig, aligned_virus, startPosition, endPosition

    def calculateCoverage(self, alignments, virusLength):
        """
        The calculateCoverage function is used to calculate this coverage. It takes as input the alignments and the length of the virus sequence. It creates a boolean array of the same length as the virus sequence, where each position in the array corresponds to a position in the virus sequence. If a position in the virus sequence is covered by an alignment, the corresponding position in the array is set to True. The coverage is then calculated as the sum of the True values in the array (i.e., the number of positions covered by alignments) divided by the length of the virus sequence, multiplied by 100 to get a percentage.
        """
        coverageArray = np.zeros(virusLength, dtype=bool)
        for alignment in alignments:
            start = alignment["startPosition"]
            end = alignment["endPosition"]
            if start != -1 and end != -1:
                coverageArray[start:end] = True
        coverage = np.sum(coverageArray) / virusLength * 100
        return coverage

    def processContigs(self):
        """
        Generator function for improved memory management
        """
        for index, contig in enumerate(self.contigs):
            for virus, virusData in self.viruses.items():
                yield index, contig, virusData["sequence"], virusData

    def processContig(self, args):
        index, contig, virusSequence, virusData = args
        maxScore, alignedContig, alignedVirus, startPosition, endPosition = (
            self.smith_waterman(contig, virusSequence)
        )
        result = {
            "contig": index + 1,
            "alignmentScore": int(maxScore),
            "alignedContigSubstring": alignedContig,
            "alignedVirusSubstring": alignedVirus,
            "lengthOfAlignment": int(len(alignedContig)),
            "startPosition": int(startPosition),
            "endPosition": int(endPosition),
        }
        return virusData["name"], result

    def search(self):
        logging.info("Search for Viruses: ")
        virusContigObj = {}
        viruses = self.viruses
        vIndex = 0

        totalContigs = len(self.contigs)
        processedContigs = 0

        for virus, virusData in viruses.items():
            print(f"Virus to align: {virus}")
            print(f"Length of virus: {len(virusData['sequence'])}")

            vStart = time.time()
            virusSequence = virusData["sequence"]
            vIndex += 1
            alignments = []

            with ProcessPoolExecutor() as executor:
                results = executor.map(self.processContig, self.processContigs())

                for virusName, result in results:
                    processedContigs += 1
                    if processedContigs == totalContigs // 4:
                        print(
                            f"25% of the contigs have been processed in {round(time.time() - vStart,3)} seconds, from the start."
                        )
                    elif processedContigs == totalContigs // 2:
                        print(
                            f"50% of the contigs have been processed in {round(time.time() - vStart,3)} seconds, from the start."
                        )
                    elif processedContigs == (totalContigs * 3) // 4:
                        print(
                            f"75% of the contigs have been processed in {round(time.time() - vStart,3)} seconds, from the start."
                        )

                    if virusName not in virusContigObj:
                        virusContigObj[virusName] = [result]
                    else:
                        virusContigObj[virusName].append(result)
            vStop = time.time()
            logging.info(f"\tVirus {virusData['name']} completed in {vStop-vStart}")
            coverage = self.calculateCoverage(alignments, len(virusSequence))
            virusContigObj[virusData["name"]].append({"coverage": coverage})
            logging.info(f"\tVirus coverage: {coverage}")

        # print(virusContigObj)
        virusContigObjFile = os.path.join(self.outputDataDir, "VirusContigObject.json")
        with open(virusContigObjFile, "w") as file:
            json.dump(virusContigObj, file)
