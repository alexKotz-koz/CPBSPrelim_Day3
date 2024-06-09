import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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
            virusAbundance[vName] = {"abundance": vProportion * 100}
        return virusAbundance

    def generateReport(self):
        virusAbundance = self.virusAbundance()
        print(virusAbundance)
        # Convert the dictionary to a pandas DataFrame
        df = pd.DataFrame(virusAbundance).T.reset_index()
        df.columns = ["Virus", "Abundance"]

        # Create the bar plot
        sns.barplot(x="Virus", y="Abundance", data=df)
        plt.show()
