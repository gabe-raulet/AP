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

    @var self.seqs: Tuple of raw sequencing reads.

    """
    def __init__(self, seqs : tuple):
        self.n = len(seqs)
        self.g = {u : dict() for u in range(self.n)}
        self.seqs = seqs

    def get_pruned(self):
        contained = set()
        for u in self.g:
            for v in self.g[u]:
                if self.g[u][v] == -1:
                    contained.add(v)

        pruned = OverlapGraph(self.seqs)

        for u in self.g:
            if not u in contained:
                for v in self.g[u]:
                    if not v in contained:
                        pruned.add_overlap(u, v, self.g[u][v])

        return pruned

    def get_meyer_tr(self, fuzz):

        VACANT = 0
        INPLAY = 1
        ELIMINATED = 2

        n = self.n
        g = self.g
        mark = [VACANT for _ in range(n)]
        reduce = {}

        for v in g:
            reduce[v] = {}
            for w in g[v]:
                reduce[v][w] = False

        neighbors = lambda v: sorted([(w, g[v][w]) for w in g[v]], key=lambda t: t[1])

        for v in range(n):
            for w in g[v]:
                mark[w] = INPLAY
            longest = max([g[v][w] for w in g[v]] + [0]) + fuzz

            v_neighbors = neighbors(v)

            for w, lvw in v_neighbors:
                if mark[w] == INPLAY:
                    w_neighbors = neighbors(w)
                    for x, lwx in w_neighbors:
                        if lvw + lwx <= longest:
                            if mark[x] == INPLAY:
                                mark[x] = ELIMINATED

            for w, lvw in v_neighbors:
                w_neighbors = neighbors(w)
                first = True
                for x, lwx in w_neighbors:
                    if lwx < fuzz or first:
                        if mark[x] == INPLAY:
                            mark[x] = ELIMINATED
                    first = False

            for w, lvw in v_neighbors:
                if mark[w] == ELIMINATED:
                    reduce[v][w] = True
                mark[w] = VACANT

        tr = OverlapGraph(self.seqs)

        for u in g:
            for w in g[u]:
                if not reduce[u][w]:
                    tr.add_overlap(u, w, g[u][w])

        return tr

    def get_naive_tr(self, fuzz):

        g = self.g
        reduce = {}

        for u in g:
            reduce[u] = {}
            for v in g[u]:
                reduce[u][v] = False

        for u in g:
            for v in g[u]:
                for w in g[v]:
                    if w in g[u]:
                        if g[u][v] + g[v][w] >= g[u][w] + fuzz:
                            reduce[u][w] = True

        tr = OverlapGraph(self.seqs)

        for u in g:
            for w in g[u]:
                if not reduce[u][w]:
                    tr.add_overlap(u, w, g[u][w])

        return tr

    def add_overlap(self, u : int, v : int, l : int):
        """
        @instance_method: add_overlap

            Adds an overlap edge to the graph.

        @param u: The source read index.

        @param v: The target read index.

        @param l: The length of the sequence that labels this overlap. This is
                     the suffix of v that overhangs past the last matching
                     symbol of u and v read from left to right. If v is
                     contained in u, then l == -1.

        @raise IndexError: If given invalid read indices.

        @raise ValueError: If the labelled suffix is the empty string.
        """

        if u < 0 or u >= self.n or v < 0 or v >= self.n:
            raise IndexError("Trying to add overlap with out of bounds indices")

        if l == 0:
            raise ValueError("Overlap suffix must not be the empty string")

        if v in self.g[u] and self.g[u][v] <= l:
            return

        self.g[u][v] = l

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
        G.vs['seq'] = self.seqs
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

        @raise ValueError: If len(seqs) != len(records).

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

        g = cls(seqs)

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
                    """
                    If vpos >= upos + ulen, then the closest read that starts
                    to the right of u is out of range, so we can cutoff early.
                    """
                    break

                if vpos + vlen <= upos + ulen:
                    """
                    If vpos + vlen <= upos + ulen, then v must be contained inside
                    u.
                    """
                    g.add_overlap(u, v, -1)
                else:
                    """
                    If we reach here, then v extends to the right of u by >= 1 nucleotides.
                    If upos == vpos, then u must be contained within v.
                    """
                    if upos == vpos:
                        g.add_overlap(v, u, -1)
                    else:
                        """
                        Reaching here means (upos < vpos) and (upos + ulen < vpos + vlen)
                        """
                        suflen = vpos + vlen - upos - ulen
                        g.add_overlap(u, v, suflen)

        return g
