import pandas as pd
import numpy as np
import logging
from logging.handlers import RotatingFileHandler

import queue
import argparse
import time
import os
import sys
import cProfile
import json

from components.importBioSample import ImportBioSample
from components.importVirus import ImportVirus
from components.qc import QualityControl

from components.readsToKmers import ReadsToKmers
from components.deBruijnGraph import DeBruijnGraph
from components.createContigs import CreateContigs

# from components.searchForViruses_SW import SearchForViruses
# from components.searchForViruses_SW_PP import SearchForViruses
from components.searchForViruses import SearchString

from components.viromeReport import ViromeReport
from components.codeReport import CodeReport

logDir = "data/logs"
os.makedirs(logDir, exist_ok=True)
# Set up logging with a rotating file handler for the main.py file
log_file = os.path.join(logDir, "app.log")

handler = RotatingFileHandler(
    log_file,
    maxBytes=1024 * 1024,  # 1MB
    backupCount=5,  # Keep 5 backup logs
)

formatter = logging.Formatter("%(message)s")
handler.setFormatter(formatter)

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(handler)


def main():
    logging.info("Main: ")
    scriptDir = os.path.dirname(__file__)
    dataDir = "data"
    dataDir = os.path.join(scriptDir, "data")
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

    componentRunTimes = {}

    viruses = {}
    virusDataDir = "./data/virus_data"
    virusFile = "sequences_20240607_3345067.fasta"
    virusFile2 = "sequences_20240607_570283.fasta"
    virusFile3 = "sequences_20240607_9774926.fasta"
    virusFile4 = "sequences_20240607_5959983.fasta"
    ncbi4VirusesFile = "ncbi_random_4_viruses.fasta"
    NCLDVFile = "sequences_ Nucleocytoviricota.fasta"
    PolBGeneFile = "kegg_polB.fasta"
    A32GeneFile = "ncbi_A32.fna"
    D5GeneFile = "ncbi_D5.fna"
    RNAplGeneFile = "kegg_RNApl.fasta"
    RNApsGeneFile = "kegg_RNAps.fasta"
    mRNAcGeneFile = "kegg_mRNAc.fasta"
    RNRSFIIGeneFile = "kegg_RNRSFII.fasta"
    VLTF3GeneFile = "kegg_VLTF3.fasta"

    # list of file locations, can be modified by adjusting
    NCLDVFileLocation = [os.path.join(virusDataDir, NCLDVFile)]
    ncbi4VirusesFileLocation = [os.path.join(virusDataDir, ncbi4VirusesFile)]
    NCLDVGeneFileLocations = [
        os.path.join(virusDataDir, PolBGeneFile),
        os.path.join(virusDataDir, A32GeneFile),
        os.path.join(virusDataDir, D5GeneFile),
        os.path.join(virusDataDir, RNAplGeneFile),
        os.path.join(virusDataDir, RNApsGeneFile),
        os.path.join(virusDataDir, mRNAcGeneFile),
        os.path.join(virusDataDir, RNRSFIIGeneFile),
        os.path.join(virusDataDir, VLTF3GeneFile),
    ]
    allVirusDataFileLocations = [
        os.path.join(virusDataDir, virusFile),
        os.path.join(virusDataDir, virusFile2),
        os.path.join(virusDataDir, virusFile3),
        os.path.join(virusDataDir, virusFile4),
        os.path.join(virusDataDir, NCLDVFile),
    ]
    syntheticVirusFileLocation = [
        os.path.join(os.path.join(dataDir, "synthetic_data"), "synthetic_virus.fasta")
    ]

    if "synthetic" in biosampleFile:
        virusDataFileLocations = syntheticVirusFileLocation
    else:
        virusDataFileLocations = ncbi4VirusesFileLocation

    if biosampleFile == "synthetic":
        biosampleFile = "synthetic_biosample.fastq"
    bioStart = time.time()
    importBioSampleInstance = ImportBioSample(biosampleFile=biosampleFile)
    biosample = importBioSampleInstance.importBioSample()
    bioStop = time.time()
    bioTotal = bioStop - bioStart
    componentRunTimes["importBioSample"] = bioTotal

    virusStart = time.time()
    importVirusInstance = ImportVirus()
    viruses = importVirusInstance.importVirusData(fileLocations=virusDataFileLocations)
    virusStop = time.time()
    virusTotal = virusStop - virusStart
    componentRunTimes["importVirus"] = virusTotal

    qcStart = time.time()
    qualityControlInstance = QualityControl(biosample=biosample)
    cleanedBiosample, minimumReadLength, qualityControlReport, qcMetadata = (
        qualityControlInstance.qualityControl()
    )
    qcStop = time.time()
    qcTotal = qcStop - qcStart
    componentRunTimes["qc"] = qcTotal

    if k > (minimumReadLength - 2):
        print(
            f"\nK must be at least one less than the size of the smallest read.\nMinimum read length for this sample is: {minimumReadLength}"
        )
        sys.exit(1)

    rtkStart = time.time()
    readsToKmersInstance = ReadsToKmers(readsData=cleanedBiosample, k=k)
    kmerPool = readsToKmersInstance.extractKmers()
    rtkStop = time.time()
    rtkTotal = rtkStop - rtkStart
    componentRunTimes["readsToKmers"] = rtkTotal
    logging.info(f"Time Stamp: Reads to Kmers finished in {rtkTotal}")
    print(f"Time Stamp: Reads to Kmers finished in {rtkTotal}")

    # Do not delete
    with open("data/logs/r-kmerPool.json", "w") as file:
        json.dump(kmerPool, file)

    dbgStart = time.time()
    debruijnGraphInstance = DeBruijnGraph(kmerPool=kmerPool, k=k)
    nodes, edges = debruijnGraphInstance.constructGraph()
    dbgStop = time.time()
    dbgTotal = dbgStop - dbgStart
    componentRunTimes["deBruijnGraph"] = dbgTotal
    logging.info(f"Time Stamp: DeBruijn Graph finished in {dbgTotal}")
    print(f"Time Stamp: DeBruijn Graph finished in {dbgTotal}")

    ccStart = time.time()
    createContigsInstance = CreateContigs(graph=edges)
    contigs = createContigsInstance.createContigs()
    ccStop = time.time()
    ccTotal = ccStop - ccStart
    componentRunTimes["createContigs"] = ccTotal
    logging.info(f"Time Stamp: Create Contigs finished in {ccTotal}")
    print(f"Time Stamp: Create Contigs finished in {ccTotal}")

    # Add comments about switch to the SW instances
    sfvStart = time.time()
    searchForVirusesInstance = SearchString(
        viruses, "data/logs/r-kmerPool.json", contigs, k
    )
    virusesInBiosample = searchForVirusesInstance.searchString()
    sfvStop = time.time()
    sfvTotal = sfvStop - sfvStart
    componentRunTimes["searchForViruses"] = sfvTotal
    logging.info(f"Time Stamp: Find Viruses finished in {sfvTotal}")
    print(f"Time Stamp: Find Viruses finished in {sfvTotal}")

    viromeReportInstance = ViromeReport(
        contigs, virusesInBiosample, biosampleFile, qcMetadata
    )
    viromeReportInstance.generateReport()

    codeReportInstance = CodeReport(
        qualityControlReport=qualityControlReport,
        contigs=contigs,
        componentRunTimes=componentRunTimes,
        qcMetadata=qcMetadata,
        biosampleFile=biosampleFile,
        k=k,
    )
    codeReportInstance.generateReport()


if __name__ == "__main__":
    start = time.time()
    main()
    stop = time.time()
    logging.info(f"Project execution time: {stop-start}")
