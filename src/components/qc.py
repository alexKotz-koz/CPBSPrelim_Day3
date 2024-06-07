import statistics
import json
import os
import pandas as pd
import logging


class QualityControl:
    def __init__(self, biosample):
        self.biosample = biosample
        self.dataDir = "./data/logs"

    # Phred and P Error calculations: https://www.drive5.com/usearch/manual/quality_score.html

    # TODO: Plot distribution of mean Q per read (ignoring ^)

    def asciiToPhred(self, character):
        return ord(character) - 33

    def calculatePError(self, q):
        return round(10 ** (-q / 10.0), 5)

    def qualityControl(self):
        biosample = self.biosample
        cleanedBiosample = {}
        qualityControlReport = {}
        qscores = []
        perrorscores = []
        sequenceLengths = []
        for id, read in biosample.items():
            id = id
            sequence = read["sequence"]
            quality = read["quality"]

            qscores = [self.asciiToPhred(q) for q in quality]

            # simulate "bases were trimmed from the end of reads if the quality score was < 20" -https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3728647/#SD2 (Supplemental Information section 5.2)
            trimmedScores = []
            trimmedSequence = []

            ### Fill in this section with the logic: Take the second half of the qscores list, check to see if the q value for each q in the second half of the qscores list is less than 20, if it is less than 20, remove the base that corresponds to that index from the sequence list
            half_index = len(qscores) // 2
            for index in range(half_index, len(qscores)):
                if qscores[index] >= 20:
                    trimmedScores.append(qscores[index])
                    trimmedSequence.append(sequence[index])

            qscores = qscores[: len(qscores) // 2] + trimmedScores
            sequence = sequence[: len(sequence) // 2] + "".join(trimmedSequence)

            sequenceLengths.append(len(sequence))

            readMedianQ = statistics.median(qscores)
            readMeanQ = statistics.mean(qscores)
            if readMedianQ >= 20:
                sequence = "".join(sequence)
                cleanedBiosample[id] = {"sequence": sequence}
                qscoreLen = len(qscores)
                qscores = []

            qualityControlReport[id] = {
                "meanQ": readMeanQ,
                "medianQ": readMedianQ,
                "length": len(sequence),
            }
        avgSeqLength = sum(sequenceLengths) / len(sequenceLengths)

        qcReportFile = os.path.join(self.dataDir, "QualityControlReport.json")
        with open(qcReportFile, "w") as file:
            json.dump(qualityControlReport, file)

        cleanedBiosampleFile = os.path.join(self.dataDir, "CleanedBioSample.json")
        with open(cleanedBiosampleFile, "w") as file:
            json.dump(cleanedBiosample, file)

        biosampleList = [{"id": id, **data} for id, data in cleanedBiosample.items()]
        biosampleDf = pd.DataFrame(biosampleList)

        logging.info("Quality Control: ")
        logging.info(f"\t# of reads: {len(cleanedBiosample)}")
        logging.info(f"\tAverage read length: {avgSeqLength}")
        logging.info(f"\tMinimum read length: {min(sequenceLengths)}")
        logging.info(f"\tMaximum read length: {max(sequenceLengths)}")

        return biosampleDf, min(sequenceLengths)
