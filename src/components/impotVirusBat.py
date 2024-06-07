import os
import logging


class ImportVirusBat:
    def __init__(self):
        self.dataDir = "./data/virus_data/bat"
        self.virusData = []

    def importVirusData(self):
        for file in os.listdir(self.dataDir):
            fileLocation = os.path.join(self.dataDir, file)
            with open(fileLocation, "r") as virusDataFile:
                self.virusData.extend(virusDataFile.readlines())

        sequence = []
        virusDataDict = {}
        virusNames = []
        for line in self.virusData:
            if line[0] == ">":

                meta = line.split(" ")
                accession = meta[0]
                name = meta[1:]
                if sequence:  # save the previous record
                    sequence = "".join(sequence)
                    virusDataDict[accession] = {
                        "name": name,
                        "sequence": sequence,
                    }
                    sequence = []

            elif line.strip():  # non-empty line
                sequence.append(line.strip())
        if sequence:  # save the last record
            sequence = "".join(sequence)
            virusDataDict[accession] = {
                "name": name,
                "sequence": sequence,
            }
        print(len(virusDataDict))

        logging.info(f"Import Viruses: ")
        logging.info(f"\tThere are {len(virusDataDict)} viruses in this dataset")
        return virusDataDict
