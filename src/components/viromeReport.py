import os


class ViromeReport:
    def __init__(self, contigs, virusesInBiosample):
        self.contigs = contigs
        self.virusesInBiosample = virusesInBiosample

    def generateReport(self):
        contigs = self.contigs
        vInB = self.virusesInBiosample
        print(len(contigs))
        # print(vInB)
