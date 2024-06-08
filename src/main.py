import pandas as pd
import numpy as np
import logging
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
from components.bacteriaRemoval import BacteriaRemoval

from components.readsToKmers import ReadsToKmers
from components.deBruijnGraph import DeBruijnGraph
from components.createContigs import CreateContigs

from components.findVirusesInReads import FindViruses
from components.searchForViruses_SW import SearchForViruses
from components.searchForViruses_SS import SearchString

# from components.searchForViruses_old import SearchForViruses


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
    importBioSampleInstance = ImportBioSample(biosampleFile=biosampleFile)
    biosample, biosampleDf = importBioSampleInstance.importBioSample()

    importVirusInstance = ImportVirus()
    viruses = importVirusInstance.importVirusData()

    qualityCheckInstance = QualityControl(biosample=biosample)
    cleanedBiosample, minimumReadLength = qualityCheckInstance.qualityControl()

    if k > (minimumReadLength - 2):
        print(
            f"\nK must be at least one less than the size of the smallest read.\nMinimum read length for this sample is: {minimumReadLength}"
        )
        sys.exit(1)

    bacteriaRemovalInstance = BacteriaRemoval()
    bacteriaRemovalInstance.removeBacteria()

    rtkStart = time.time()
    readsToKmersInstance = ReadsToKmers(readsData=cleanedBiosample, k=k)
    kmerPool, _ = readsToKmersInstance.extractKmers()
    rtkStop = time.time()
    logging.info(f"Time Stamp: Reads to Kmers finished in {rtkStop-rtkStart}")
    print(f"Time Stamp: Reads to Kmers finished in {rtkStop-rtkStart}")

    dbgStart = time.time()
    debruijnGraphInstance = DeBruijnGraph(kmerPool=kmerPool, k=k)
    nodes, edges, orphanedNodes, orphanedSubgraphs = (
        debruijnGraphInstance.constructGraph()
    )
    print(f"Orphaned Nodes: {orphanedNodes}, Orphaned Subgraphs: {orphanedSubgraphs}")
    dbgStop = time.time()
    logging.info(f"Time Stamp: DeBruijn Graph finished in {dbgStop-dbgStart}")
    print(f"Time Stamp: DeBruijn Graph finished in {dbgStop-dbgStart}")

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

    """fvStart = time.time()
    findVirusesInstance = FindViruses(reads=cleanedBiosample, viruses=viruses)
    findVirusesInstance.findViruses()
    fvStop = time.time()
    logging.info(f"Time Stamp: Find Viruses finished in {fvStop-fvStart}")
    print(f"Time Stamp: Find Viruses finished in {fvStop-fvStart}")"""

    sfvStart = time.time()
    searchForVirusesInstance = SearchForViruses(viruses=viruses, contigs=contigs, k=k)
    searchForVirusesInstance.search()
    sfvStop = time.time()
    logging.info(f"Time Stamp: Find Viruses finished in {sfvStop-sfvStart}")
    print(f"Time Stamp: Find Viruses finished in {sfvStop-sfvStart}")

    """sfvStart = time.time()
    for id, virus in viruses.items():
        searchStringInstance = SearchString(
            virus=virus, readsKmerPool=kmerPool, contigs=contigs, k=k
        )
        searchStringInstance.searchString()
        break
    sfvStop = time.time()
    logging.info(f"Time Stamp: Find Viruses finished in {sfvStop-sfvStart}")
    print(f"Time Stamp: Find Viruses finished in {sfvStop-sfvStart}")"""


if __name__ == "__main__":
    main()
