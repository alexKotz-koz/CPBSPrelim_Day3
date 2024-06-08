import numpy as np
import os
import json
import logging
import time
from utils.utils import ToKmers


class FindViruses:
    def __init__(self, reads, viruses):
        self.reads = reads
        self.viruses = viruses
        self.dataDir = "./data/output_data"

    #'>OK295287.1': {'name': 'Severe acute respiratory syndrome coronavirus 2', 'sequence':

    # 1. Initialization
    # 2. Matrix Filling
    # 3. Trace Back

    def smith_waterman(
        self, seq1, seq2, match_score=3, mismatch_score=-1, gap_penalty=-2
    ):
        """
        Perform Smith-Waterman alignment on two sequences.

        Args:
        - seq1 (str): First sequence to align.
        - seq2 (str): Second sequence to align.
        - match_score (int): Score for a match.
        - mismatch_score (int): Score for a mismatch.
        - gap_penalty (int): Penalty for opening a gap.

        Returns:
        - alignment_score (int): Score of the best alignment.
        - aligned_seq1 (str): First sequence with gaps for alignment.
        - aligned_seq2 (str): Second sequence with gaps for alignment.
        """
        # Initialize scoring matrix
        rows = len(seq1) + 1
        cols = len(seq2) + 1
        score_matrix = np.zeros((rows, cols), dtype=int)

        # Initialize traceback matrix
        traceback_matrix = np.zeros((rows, cols), dtype=int)

        # Fill in scoring and traceback matrices
        for i in range(1, rows):
            for j in range(1, cols):
                match = score_matrix[i - 1][j - 1] + (
                    match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
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
        aligned_seq1 = ""
        aligned_seq2 = ""
        i, j = max_index
        while i > 0 and j > 0 and score_matrix[i][j] != 0:
            if traceback_matrix[i][j] == 1:
                aligned_seq1 = seq1[i - 1] + aligned_seq1
                aligned_seq2 = seq2[j - 1] + aligned_seq2
                i -= 1
                j -= 1
            elif traceback_matrix[i][j] == 2:
                aligned_seq1 = seq1[i - 1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                i -= 1
            else:
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[j - 1] + aligned_seq2
                j -= 1
        print(f"Done: {max_score} {aligned_seq1} {aligned_seq2}")
        return max_score, aligned_seq1, aligned_seq2

    def hammingDistance(self, virus, contig):
        distance = 0
        for index in range(len(virus)):
            if contig[index] != virus[index]:
                distance += 1
        return distance

    def contigInVirus(self, virus, contig):
        if contig in virus:
            print("yes")

    def findViruses(self):
        contigsInVirus = {}
        viruses = self.viruses
        reads = self.reads
        matchLenContigVirus = {}
        logging.info("Find Viruses in Reads: ")
        for key, virus in viruses.items():
            vName = virus["name"]
            vSeq = virus["sequence"]
            sameLen = 0
            for index, read in reads.iterrows():
                id = read["id"]
                rSeq = read["sequence"]
                if len(rSeq) == len(vSeq):
                    sameLen += 1
                    distance = self.hammingDistance(vSeq, rSeq)
                    coverage = (len(rSeq) - distance) / len(rSeq)

                    if key not in matchLenContigVirus:
                        matchLenContigVirus[vName] = [
                            {
                                "read": index,
                                "distance": distance,
                                "sequenceLen": len(rSeq),
                                "coverage": coverage,
                            }
                        ]
                    else:
                        matchLenContigVirus[vName].append(
                            {
                                "read": index,
                                "distance": distance,
                                "sequenceLen": len(rSeq),
                                "coverage": coverage,
                            }
                        )
            for virus, data in matchLenContigVirus.items():
                for read in data:
                    if read["distance"] < 100:
                        print(virus)
        # logging.info(f"\t{matchLenContigVirus}")
        #    start = time.time()
        #    stop = time.time()
        #    # print(f"{virus['name']} done in {stop-start}")

        contigsInVirusFile = os.path.join(self.dataDir, "ContigsInVirus.json")
        with open(contigsInVirusFile, "w") as file:
            json.dump(matchLenContigVirus, file)
