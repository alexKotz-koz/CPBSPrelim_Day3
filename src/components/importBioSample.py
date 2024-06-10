import pandas as pd
import os


class ImportBioSample:
    def __init__(self, biosampleFile):
        self.biosampleFile = biosampleFile
        self.biosample = {}

    def importBioSample(self):
        dataDir = "./data/biosample_data"
        os.makedirs(dataDir, exist_ok=True)
        biosample = self.biosample
        biosampleFile = self.biosampleFile
        biosampleDataDir = ""

        if "biofilm" in biosampleFile:
            biosampleDataDir = os.path.join(dataDir, "biofilm")
        elif "cryoconite" in biosampleFile:
            biosampleDataDir = os.path.join(dataDir, "cryoconite")
        elif "bat" in biosampleFile:
            biosampleDataDir = os.path.join(dataDir, "bat")
        elif "synthetic" in biosampleFile:
            biosampleDataDir = dataDir
        elif "test" in biosampleFile:
            biosampleDataDir = ""

        if biosampleDataDir == "":
            biosampleFileLocation = "test_data/testbiosample.fastq"
        else:
            biosampleFileLocation = os.path.join(biosampleDataDir, biosampleFile)

        with open(biosampleFileLocation, "r") as biosampleFile:
            while True:
                idline = biosampleFile.readline().split(" ")
                if not idline[0]:
                    break
                id = idline[0]
                seq = biosampleFile.readline().strip()
                plus = biosampleFile.readline().strip()
                quality = biosampleFile.readline().strip()

                biosample[id] = {"sequence": seq, "quality": quality}

        biosampleList = [{"id": id, **data} for id, data in biosample.items()]

        return biosample
