import os
import logging


class ImportVirus:
    def __init__(self):
        self.dataDir = "./data/virus_data"
        self.virusData = []

    def importVirusData(self, fileLocations):
        for fileLocation in fileLocations:
            with open(fileLocation, "r") as virusDataFile:
                self.virusData.extend(virusDataFile.readlines())
        """with open(os.path.join(self.dataDir, self.NCLDVFile), "r") as file:
            self.virusData.extend(file.readlines())"""
        sequence = []
        virusDataDict = {}
        virusNames = []
        for line in self.virusData:
            if line[0] == ">":

                metadata = line.split("|")
                accession = metadata[0].strip()
                organismName = metadata[1].strip()

                if organismName in virusNames:
                    continue

                virusNames.append(organismName)

                if sequence:  # save the previous record
                    sequence = "".join(sequence)
                    virusDataDict[accession] = {
                        "name": organismName,
                        "sequence": sequence,
                    }
                    sequence = []

            elif line.strip():  # non-empty line
                sequence.append(line.strip())
        if sequence:  # save the last record
            sequence = "".join(sequence)
            virusDataDict[accession] = {
                "name": organismName,
                "sequence": sequence,
            }
        virusDataDict = dict(
            sorted(virusDataDict.items(), key=lambda item: len(item[1]["sequence"]))
        )
        logging.info(f"Import Viruses: ")
        logging.info(f"\tThere are {len(virusDataDict)} viruses in this dataset")
        return virusDataDict
