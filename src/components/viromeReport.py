import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class ViromeReport:
    def __init__(self, contigs, virusesInBiosample):
        self.contigs = contigs
        self.virusesInBiosample = virusesInBiosample
        self.reportDir = "data/reports"
        os.makedirs(self.reportDir, exist_ok=True)

    def virusAbundance(self):
        virusAbundance = {}
        contigs = self.contigs
        totalNumContigs = len(contigs)
        vInB = self.virusesInBiosample
        for virus in vInB:
            vName = virus["virus"]
            numContigs = virus["numContigsInVirus"]
            vProportion = numContigs / totalNumContigs
            virusAbundance[vName] = {"abundance": vProportion * 100}
        return virusAbundance

    def generateReport(self):
        # virus abundance
        virusAbundance = self.virusAbundance()
        df = pd.DataFrame(virusAbundance).T.reset_index()
        df.columns = ["Virus", "Abundance"]
        sns.barplot(x="Virus", y="Abundance", data=df)
        plt.show()
