import pandas as pd
import os


class ImportBioSample:
    def __init__(self, biosampleFile):
        self.biosampleFile = biosampleFile
        self.biosample = {}

    def importBioSample(self):
        scriptDir = os.path.dirname(os.path.dirname(__file__))
        dataDir = "data"
        dataDir = os.path.join(scriptDir, "data")
        self.logsDataDir = os.path.join(dataDir, "logs")
        self.outputDataDir = os.path.join(dataDir, "output_data")
        biosample_dataDir = os.path.join(dataDir, "biosample_data")
        biosample = self.biosample
        biosampleFile = self.biosampleFile
        biosampleDataDir = ""

        if "biofilm" in biosampleFile:
            biosampleDataDir = os.path.join(biosample_dataDir, "biofilm")
        elif "cryoconite" in biosampleFile:
            biosampleDataDir = os.path.join(biosample_dataDir, "cryoconite")
        elif "synthetic" in biosampleFile:
            biosampleDataDir = os.path.join(dataDir, "synthetic_data")
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
