from igraph import Graph as IGraph

class OverlapGraph(object):
    """
    @class OverlapGraph:

        Class for an overlap graph G = (V,E) where V = \{0,...,self.n-1}
        is the set of reads/vertices and E is the set of read overlaps.

    @var self.n: Number of reads/vertices in graph.

    @var self.g: Adjacency list of graph. For a given u \in V,
                 if u -> v  is an edge in the graph then self.g[u][v]
                 exists and is equal to the labelled suffix of the
                 read overlap u -> v.
    """
    def __init__(self, n : int):
        self.n = n
        self.g = {u : dict() for u in range(n)}

    def add_overlap(self, u : int, v : int, vsuf : str):
        """
        @instance_method: add_overlap

            Adds an overlap edge to the graph.

        @param u: The source read index.

        @param v: The target read index.

        @param vsuf: The sequence that labels this overlap. This is
                     the suffix of v that overhangs past the last matching
                     symbol of u and v read from left to right.

        @raise IndexError: If given invalid read indices.

        @raise ValueError: If the labelled suffix is the empty string.
        """

        if u < 0 or u >= self.n or v < 0 or v >= self.n:
            raise IndexError("Trying to add overlap with out of bounds indices")

        if len(vsuf) == 0:
            raise ValueError("Overlap suffix must not be the empty string")

        if v in self.g[u] and len(self.g[u][v]) <= len(vsuf):
            return

        self.g[u][v] = vsuf


    def num_edges(self) -> int:
        """
        @instance_method num_edges:

            Calculates the number of overlaps/edges in the graph.

        This operation takes O(|V|) time each call because the edge count
        is not stored internally and must be recalculated each time.
        """
        return sum(len(adj) for v,adj in self.g.items())

    def get_igraph(self) -> IGraph:
        """
        @instance_method get_igraph:

            Constructs and returns a new igraph object created from the contents
            of the self graph.

        @returns G: igraph object
        """
        G = IGraph(self.n, directed = True)
        edges = []
        sufs = []
        for u, adj in self.g.items():
            for v, suffix in adj.items():
                edges.append((u,v))
                sufs.append(suffix)
        G.add_edges(edges)
        G.es['suf'] = sufs
        return G

    @classmethod
    def generate(cls, seqs : tuple, records : list, genome_length : int): # -> cls: (annotations don't work for this use case apparently)
        """
        @classmethod generate:

            Constructs a new OverlapGraph object usingt input reads and
            their gold standard mapping back to original reference.

        @param seqs: Tuple of reads ordered be read index.

        @param records: Read triplets. Look at AssemblyProblem.create_reads.__doc__ function for definition.

        @param genome_length: Reference genome length; used to calculate overlaps
                              that cross the circular barrier.

        @raise TypeError: If input parameters have wrong type.

        @raise ValueErorr: If len(seqs) != len(records).

        @raise AssertionError: If records are incorrectly sorted (should never happen).

        @return instantiated OverlapGraph object.
        """

        if not isinstance(seqs, tuple):
            raise TypeError("Input sequences are not a tuple.")

        if not isinstance(records, list):
            raise TypeError("Input records not a list.")

        if not isinstance(genome_length, int):
            raise TypeError("Genome length not an integer.")

        n = len(seqs)

        if len(records) != n:
            raise ValueError("len(seqs) != len(records)")

        g = cls(n)

        records = sorted(records, key=lambda t: t[1]) * 2

        for ufind in range(n):

            u, upos, urev = records[ufind]
            useq = seqs[u]
            ulen = len(useq)

            for vfind in range(ufind+1, 2*n):

                if vfind < n:
                    v, vpos, vrev = records[vfind]
                else:
                    v, vpos, vrev = records[vfind - n]
                    vpos += genome_length

                assert upos <= vpos

                vseq = seqs[v]
                vlen = len(vseq)

                if vpos >= upos + ulen:
                    break

                if vpos + vlen > upos + ulen:
                    suflen = vpos + vlen - upos - ulen
                    g.add_overlap(u, v, vseq[-suflen : ])

        return g


