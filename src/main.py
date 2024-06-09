import pandas as pd
import numpy as np
import logging
import logging.handlers
import queue
import argparse
import time
import os
import sys
import cProfile
import pstats

from components.importBioSample import ImportBioSample
from components.importVirus import ImportVirus
from components.impotVirusBat import ImportVirusBat
from components.qc import QualityControl

from components.readsToKmers import ReadsToKmers
from components.deBruijnGraph import DeBruijnGraph
from components.createContigs import CreateContigs

# from components.searchForViruses_SW import SearchForViruses
# from components.searchForViruses_SW_PP import SearchForViruses
from components.searchForViruses_hamming import SearchForVirusesHamming
from components.searchForViruses_old import SearchString

from components.viromeReport import ViromeReport


logDir = "data/logs"
os.makedirs(logDir, exist_ok=True)

logging.basicConfig(
    filename=f"{logDir}/app.log",
    filemode="w",
    format="%(message)s",
    level=logging.INFO,
)


def main():
    logging.info("Main: ")
    currDir = os.getcwd()
    rootDir = os.path.dirname(currDir)
    srcDir = os.path.join(rootDir, "src")
    dataDir = os.path.join(srcDir, "data")
    # arg setup and management
    parser = argparse.ArgumentParser(
        description="Metagenomic-based virome characterizer"
    )
    parser.add_argument(
        "-biosample", type=str, help="Fastq biosample file", required=True
    )
    parser.add_argument("-k", type=int, help="size of kmer", required=True)

    args = parser.parse_args()

    biosampleFile = args.biosample
    k = args.k
    logging.info(f"\tBioSample File: {biosampleFile}")
    logging.info(f"\tSize of K = {k}")
    # Simulated Data
    viruses = {}

    if biosampleFile == "synthetic":

        if os.path.exists(srcDir):
            pass
        else:
            rootDir = os.path.join(rootDir, "CPBSPrelim_Day3")
            srcDir = os.path.join(srcDir, rootDir)
            dataDir = os.path.join(srcDir, dataDir)
            if os.path.exists(os.path.join(dataDir, "synthetic_data")):
                pass
            else:
                raise Exception(f"Syntheticdir still DNE: {dataDir}")
        syntheticDataDir = os.path.join(dataDir, "synthetic_data")
        syntheticBiosampleFile = os.path.join(syntheticDataDir, "biosample.fastq")
        syntheticVirusFile = os.path.join(syntheticDataDir, "virus.fasta")

        biosample_dict = {}
        with open(syntheticBiosampleFile, "r") as file:
            data = file.readlines()
        for index, line in enumerate(data):
            if line[0] == "@":
                id = line.strip()
                if index + 1 < len(data):
                    read = data[index + 1].strip()
                    biosample_dict[id] = read

        cleanedBiosample = pd.DataFrame(
            list(biosample_dict.items()), columns=["id", "sequence"]
        )

        with open(syntheticVirusFile, "r") as file:
            vdata = file.readlines()
        for index, item in enumerate(vdata):
            if item[0] == ">":
                id = item.strip()
                seq = vdata[index + 1].strip()
            viruses[id] = {"name": id, "sequence": seq}
    else:
        virusDataDir = "./data/virus_data"
        virusFile = "sequences_20240607_3345067.fasta"
        virusFile2 = "sequences_20240607_570283.fasta"
        virusFile3 = "sequences_20240607_9774926.fasta"
        virusFile4 = "sequences_20240607_5959983.fasta"
        NCLDVFile = "sequences_ Nucleocytoviricota.fasta"
        NCLDVFileLocation = [os.path.join(virusDataDir, NCLDVFile)]
        virusDataFileLocations = [
            os.path.join(virusDataDir, virusFile),
            os.path.join(virusDataDir, virusFile2),
            os.path.join(virusDataDir, virusFile3),
            os.path.join(virusDataDir, virusFile4),
            os.path.join(virusDataDir, NCLDVFile),
        ]

        importBioSampleInstance = ImportBioSample(biosampleFile=biosampleFile)
        biosample, biosampleDf = importBioSampleInstance.importBioSample()

        importVirusInstance = ImportVirus()
        viruses = importVirusInstance.importVirusData(fileLocations=NCLDVFileLocation)

        qualityCheckInstance = QualityControl(biosample=biosample)
        cleanedBiosample, minimumReadLength = qualityCheckInstance.qualityControl()
        if k > (minimumReadLength - 2):
            print(
                f"\nK must be at least one less than the size of the smallest read.\nMinimum read length for this sample is: {minimumReadLength}"
            )
            sys.exit(1)

    rtkStart = time.time()
    readsToKmersInstance = ReadsToKmers(readsData=cleanedBiosample, k=k)
    kmerPool, _ = readsToKmersInstance.extractKmers()
    rtkStop = time.time()
    logging.info(f"Time Stamp: Reads to Kmers finished in {rtkStop-rtkStart}")
    print(f"Time Stamp: Reads to Kmers finished in {rtkStop-rtkStart}")

    dbgStart = time.time()
    debruijnGraphInstance = DeBruijnGraph(kmerPool=kmerPool, k=k)
    nodes, edges = debruijnGraphInstance.constructGraph()
    dbgStop = time.time()
    logging.info(f"Time Stamp: DeBruijn Graph finished in {dbgStop-dbgStart}")
    print(f"Time Stamp: DeBruijn Graph finished in {dbgStop-dbgStart}")
    # print(f"edges: {edges}")
    ccStart = time.time()
    createContigsInstance = CreateContigs(graph=edges)
    # contigs, allPaths = cProfile.run(createContigsInstance.createContigs(), "output.dat")
    contigs = createContigsInstance.createContigs()
    ccStop = time.time()
    logging.info(f"Time Stamp: Create Contigs finished in {ccStop-ccStart}")
    print(f"Time Stamp: Create Contigs finished in {ccStop-ccStart}")

    """p = pstats.Stats("output.dat")
    p.sort_stats("cumulative").print_stats(
        10
    )  # Print the 10 most time-consuming functions"""

    searchForVirusesInstance = SearchString(viruses, kmerPool, contigs, k)
    virusesInBiosample = searchForVirusesInstance.searchString()

    """sfvStart = time.time()
    searchForVirusesInstance = SearchForViruses(
        viruses=viruses, contigs=cleanedBiosample, k=k
    )
    searchForVirusesInstance.search()
    sfvStop = time.time()
    logging.info(f"Time Stamp: Find Viruses finished in {sfvStop-sfvStart}")
    print(f"Time Stamp: Find Viruses finished in {sfvStop-sfvStart}")"""

    """sfvStart = time.time()
    searchForVirusesInstance = SearchForVirusesHamming(
        viruses=viruses, reads=cleanedBiosample, k=k
    )
    searchForVirusesInstance.search()
    sfvStop = time.time()
    logging.info(f"Time Stamp: Find Viruses finished in {sfvStop-sfvStart}")
    print(f"Time Stamp: Find Viruses finished in {sfvStop-sfvStart}")"""

    viromeReportInstance = ViromeReport(contigs, virusesInBiosample)
    viromeReportInstance.generateReport()


if __name__ == "__main__":
    main()
