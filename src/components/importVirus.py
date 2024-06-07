import os
import logging


class ImportVirus:
    def __init__(self):
        virusFile = "sequences_20240607_3345067.fasta"
        virusFile2 = "sequences_20240607_570283.fasta"
        virusFile3 = "sequences_20240607_9774926.fasta"
        virusFile4 = "sequences_20240607_5959983.fasta"
        NCLDVFile = "sequences_ Nucleocytoviricota.fasta"
        dataDir = "./data/virus_data"
        self.virusDataFileLocations = [
            os.path.join(dataDir, virusFile),
            os.path.join(dataDir, virusFile2),
            os.path.join(dataDir, virusFile3),
            os.path.join(dataDir, virusFile4),
            os.path.join(dataDir, NCLDVFile),
        ]
        self.virusData = []

    def importVirusData(self):
        for fileLocation in self.virusDataFileLocations:
            with open(fileLocation, "r") as virusDataFile:
                self.virusData.extend(virusDataFile.readlines())

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
        print(len(virusNames))

        logging.info(f"Import Viruses: ")
        logging.info(f"\tThere are {len(virusDataDict)} viruses in this dataset")
        return virusDataDict
