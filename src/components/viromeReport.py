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
            virusAbundance[vName] = {"abundance": vProportion}
        return virusAbundance

    def generateReport(self):
        virusAbundance = self.virusAbundance()
        print(virusAbundance)
