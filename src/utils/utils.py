import random
import string
from collections import defaultdict


def ToKmers(k, df=False, obj=False, sequence=False):

    if df and object == False:
        kmerPool = defaultdict(lambda: defaultdict(list))
        readsData = df
        for read in readsData.itertuples():
            id = read.id
            sequence = read.sequence
            kmers = [sequence[i : i + k] for i in range(len(sequence) - k + 1)]
            for index, kmer in enumerate(kmers):
                kmerPool[kmer][id].append({index: index + k})
        return kmerPool, k

    elif obj and df == False:
        kmerPool = {}
        for key, virus in obj.items():
            sequence = virus["sequence"]
        for index, base in enumerate(sequence):
            kmer = sequence[index : index + k]
            if len(kmer) >= k:
                if kmer not in kmerPool:
                    kmerPool[kmer] = [{index: index + k}]
                else:
                    kmerPool[kmer].append({index: index + k})
        return kmerPool


def generateRandomData():
    bases = ["A", "G", "C", "T", "N"]
    weights = [0.225, 0.225, 0.225, 0.225, 0.1]  # specify your weights here
    length = 300
    random.seed(2038)

    sequences = {}

    for i in range(100):
        sequence = []
        id = i + 1
        for i in range(length):
            base = random.choices(bases, weights=weights, k=1)[0]
            sequence.append(base)
        name = f"@{id}|{''.join(random.choices(string.ascii_letters, k=5))}"
        sequence = "".join(sequence)
        sequences[id] = {
            "name": name,
            "sequence": sequence,
        }

    sequences = pd.DataFrame.from_dict(
        sequences, columns=["name", "sequence"], orient="index"
    )
    sequences.to_csv("dummy_reads_100.csv", index=True, index_label="id")
