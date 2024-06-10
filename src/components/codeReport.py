import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class CodeReport:
    def __init__(self, qualityControlReport, contigs, componentRunTimes):
        self.qcReport = qualityControlReport
        self.contigs = contigs
        self.componentRunTimes = componentRunTimes
        self.reportDir = "data/reports"
        os.makedirs(self.reportDir, exist_ok=True)

    # Contigs generated
    # Size of k
    # Time for each component completion

    def qcReportAnalysis(self):
        qcReportDf = pd.DataFrame(self.qcReport).T
        fig, axs = plt.subplots(3, figsize=(10, 10))

        yLabels = [
            "Mean Q-Score (per read)",
            "Median Q-Score (per read)",
            "Length (per read)",
        ]
        colors = ["blue", "green", "red"]
        for index, column in enumerate(qcReportDf.columns):
            sns.lineplot(
                data=qcReportDf[column], ax=axs[index], marker="o", color=colors[index]
            )
            axs[index].set_xticklabels([])
            axs[index].set_ylabel(yLabels[index], weight="bold")

        fig.text(0.5, 0.04, "Reads", ha="center", va="center", weight="bold")

        fig.suptitle("Quality Control Report", weight="bold")
        plt.subplots_adjust(bottom=0.1)

        plt.show()
        return len(self.qcReport)

    def assemblyStatistics(self):
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11020464/
        # https://www.nature.com/articles/s41467-022-35735-y fig 1
        # num contigs
        # contig length distribution
        # time per component

        pass

    def generateReport(self):
        reportFile = "code_report.txt"
        numReads = self.qcReportAnalysis()
