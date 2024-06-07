import pandas as pd
import os


class ImportBioSample:
    def __init__(self, biosampleFile):
        dataDir = "./data/biosample_data"
        if "biofilm" in biosampleFile:
            biosampleDataDir = os.path.join(dataDir, "biofilm")
        elif "cryoconite" in biosampleFile:
            biosampleDataDir = os.path.join(dataDir, "cryoconite")
        elif "bat" in biosampleFile:
            biosampleDataDir = os.path.join(dataDir, "bat")
        self.biosampleFileLocation = os.path.join(biosampleDataDir, biosampleFile)

        self.biosample = {}

    def importBioSample(self):
        biosample = self.biosample
        with open(self.biosampleFileLocation, "r") as biosampleFile:
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

        biosampleDf = pd.DataFrame(biosampleList)

        return biosample, biosampleDf
