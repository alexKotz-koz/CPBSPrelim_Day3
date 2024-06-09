import os


class ViromeReport:
    def __init__(self, contigs, virusesInBiosample):
        self.contigs = contigs
        self.virusesInBiosample = virusesInBiosample

    def virusAbundance(self):
        virusAbundance = {}
        contigs = self.contigs
        totalNumContigs = len(contigs)
        vInB = self.virusesInBiosample
        for virus in vInB:
            vName = virus["virus"]
            numContigs = virus["numContigsInVirus"]
            vProportion = numContigs / totalNumContigs
            print(vProportion)

    def generateReport(self):
        self.virusAbundance()
