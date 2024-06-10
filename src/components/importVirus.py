import os
import logging


class ImportVirus:
    def __init__(self):
        self.dataDir = "./data/virus_data"
        os.makedirs(self.dataDir, exist_ok=True)
        self.virusData = []

    def importVirusData(self, fileLocations):
        # import virus files
        for fileLocation in fileLocations:
            with open(fileLocation, "r") as virusDataFile:
                self.virusData.extend(virusDataFile.readlines())

        sequence = []
        virusDataDict = {}
        virusNames = []
        for line in self.virusData:
            if line[0] == ">":
                if sequence:  # save the previous record
                    sequence = "".join(sequence).upper()
                    virusDataDict[accession] = {
                        "name": organismName,
                        "sequence": sequence,
                        "length": len(sequence),
                    }
                    sequence = []

                metadata = line.split("|")
                accession = metadata[0].strip()
                organismName = metadata[1].strip()

                # skip possible duplicate virus
                if organismName in virusNames:
                    print(f"Duplicate Virus, skipping: {organismName}")
                    continue

                virusNames.append(organismName)

            elif line.strip():  # non-empty line
                sequence.append(line.strip())

        if sequence:  # save the last record
            sequence = "".join(sequence).upper()
            virusDataDict[accession] = {
                "name": organismName,
                "sequence": sequence,
                "length": len(sequence),
            }

        virusDataDict = dict(
            sorted(virusDataDict.items(), key=lambda item: len(item[1]["sequence"]))
        )

        logging.info(f"Import Viruses: ")
        logging.info(f"\tThere are {len(virusDataDict)} viruses in this dataset")
        return virusDataDict
