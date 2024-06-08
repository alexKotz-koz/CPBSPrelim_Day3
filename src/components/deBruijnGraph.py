import pandas as pd
import numpy as np


class DeBruijnGraph:
    def __init__(self, kmerPool, k):
        self.kmerPool = kmerPool
        self.k = k

    # Input: individual kmer
    # Output: prefix and suffix of kmer
    def getPrefixSuffix(self, kmer):
        # get the length of the prefix/suffix
        length = self.k - 1
        # get the first <length> characters of the kmer for the prefix
        prefix = kmer[:length]
        # get the last <length> characters of the kmer for the suffix
        suffix = kmer[-length:]
        return prefix, suffix

    def findOrphanedNodes(self, nodes, edges):
        orphanedNodes = [
            node
            for node in nodes
            if node not in edges
            and not any(node in targets for targets in edges.values())
        ]
        return len(orphanedNodes)

    def findOrphanedSubgraphs(self, nodes, edges):
        visited = set()
        orphanedSubgraphs = 0

        for node in nodes:
            if node not in visited:
                # Start a DFS from this node
                stack = [node]
                connectedComponent = set()

                while stack:
                    currentNode = stack.pop()
                    visited.add(currentNode)
                    connectedComponent.add(currentNode)

                    if currentNode in edges:
                        for target in edges[currentNode]:
                            if target not in visited:
                                stack.append(target)

                # Check if this connected component is connected to any other components
                if not any(node in edges for node in connectedComponent) and not any(
                    node in targets
                    for targets in edges.values()
                    for node in connectedComponent
                ):
                    orphanedSubgraphs += 1

        return orphanedSubgraphs

    # Input: kmerPool (contains all unique kmers found in the reads, the read id's where each kmer exists, and the position of each kmer in each read)
    # Output: nodes (a set containing all unique nodes), edges (a list containing tuples of source->target nodes)
    def constructGraph(self):
        kmerPool = self.kmerPool
        edges = {}
        nodes = set()

        lengthOfPool = len(kmerPool)

        # get the prefix and suffix of each kmer and add to the nodes and edges data structures
        for kmer in kmerPool:
            prefix, suffix = self.getPrefixSuffix(kmer)
            if not edges.get(prefix):
                edges[prefix] = [suffix]
            else:
                edges[prefix].append(suffix)
            nodes.add(prefix)
            nodes.add(suffix)

        orphanedNodes = self.findOrphanedNodes(nodes, edges)
        orphanedSubgraphs = self.findOrphanedSubgraphs(nodes, edges)

        return nodes, edges, orphanedNodes, orphanedSubgraphs
