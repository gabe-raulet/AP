from igraph import Graph as IGraph

class OverlapGraph(object):
    """
    @class OverlapGraph:

        Class for a bidirected overlap graph G = (V,E) where V = \{0,...,self.n-1\}
        is the set of reads/vertices and E is the set of read overlaps.

        Below are examples/descriptions of each kind of overlap and their associated
        semantic representation in the graph.

                        012*45678*
        ############  u --------->         beg(u) = 3, end(u) = 9, len(u) = 9
        # Regular  #       ||||||
        # Dovetail #     v ---------->     beg(v) = 0, end(v) = 6, len(v) = 10
        ############       *12345*789

        In the above case, a suffix of u aligns to a prefix of v. This is called
        a Regular Dovetail. In the overlap graph, we represent this with a directed
        bidirected edge (u >--> v). We label this edge (semantically) with the overhanging
        suffix of v which does not align with u, namely v[end(v)..).
            Technically, this is all the information required for this scenario and this
        is a valid bidirected edge. However, since we will be using the graph to perform
        walks/traversals to build contigs, we can make the job easier right now to save
        time later. Specifically, imagine that in our traversal we passed from v to u,
        which is the opposite direction of the edge's orientation. In such a case we don't
        care about the suffix of v; we need the prefix of u, and more specifically, the
        reverse complement of the prefix of u. It is possible to calculate this from
        the information given so far. Namely, if l(u >--> v) denotes the length of the
        sequence labeling the edge (u >--> v), so that v[|v|-l..) is the label, then we
        can calculate the length of the overhanging prefix of u to be m = |u| + l - |v|,
        and then the prefix of interest is comp(u[..m)). However, it is easier to just label
        this immediately. We accomplish this by creating a shadow edge for the prefix direction,
        which can be looked up as if we were looking up the "directed edge" (v --> u), and in
        that location we store comp(u[..m)). A similar procedure of storing "shadow" edges for
        quickly looking up forward and backward traversals of a bidirected edge will be taken
        for the other kinds of edges as well.

                        012*45678*
        ############  u --------->         beg(u) = 3, end(u) = 9, len(u) = 9
        # Prefix   #       ||||||
        # Dovetail #    v <----------      beg(v) = 9, end(v) = 3, len(v) = 10
        ############       *87654*210

        In this case, a suffix of u aligns to a prefix of comp(v), which is the reverse
        complement of a suffix of v. In the overlap graph, we represent this with an
        introverted bidirected edge (u >--< v). We label the (u->v) direction with the
        sequence comp(v[..end(v)]) and the (v->u) direction with the sequence u[..beg(u)).

        TODO: finish

    @var self.n: Number of reads/vertices in graph.

    @var self.g: Adjacency list of graph. self.g[u][v] = (dir(u,v), suflen(u,v))

    #  @var self.g: Adjacency list of graph. For a given u \in V,
                 #  if u -> v  is an edge in the graph then self.g[u][v]
                 #  exists and is equal to the

    @var self.seqs: Tuple of raw sequencing reads.

    """
    def __init__(self, seqs : tuple):
        self.n = len(seqs)
        self.g = {u : dict() for u in range(self.n)}
        self.seqs = seqs

    def get_pruned(self):
        """
        @func get_pruned: TODO
        """
        contained = set()
        for u in self.g:
            for v in self.g[u]:
                if self.g[u][v][0] == -1:
                    contained.add(v)

        pruned = OverlapGraph(self.seqs)

        for u in self.g:
            if not u in contained:
                for v in self.g[u]:
                    if not v in contained:
                        pruned.add_overlap(u, v, self.g[u][v][0], self.g[u][v][1])

        return pruned

    def get_naive_tr(self, fuzz):
        """
        @func get_naive_tr: TODO
        """

        g = self.g
        reduce = {u : {} for u in g}

        for u in g:
            for v in g[u]:
                reduce[u][v] = False

        def edir(u, v): return g[u][v][0]
        def suflen(u, v): return g[u][v][1]
        def arrows(u, v): return ((g[u][v][0]>>1)&1, g[u][v][0]&1)

        for u in g:
            for v in g[u]:
                uv_head, uv_tail = arrows(u,v)
                for w in set(g[v]).intersection(g[u]):
                    vw_head, vw_tail = arrows(v,w)
                    uw_head, uw_tail = arrows(u,w)
                    if uv_head == uw_head and vw_tail == uw_tail and uv_tail != vw_head:
                        if suflen(u,w) + suflen(w,v) >= suflen(u,w) + fuzz:
                            reduce[u][w] = True
                            reduce[w][u] = True

        tr = OverlapGraph(self.seqs)

        for u in g:
            for w in g[u]:
                if not reduce[u][w]:
                    assert not reduce[w][u]
                    tr.add_overlap(u, w, edir(u, w), suflen(u, w))
                    tr.add_overlap(w, u, edir(w, u), suflen(w, u))

        return tr

    def add_overlap(self, u : int, v : int, d : int, l : int):
        """
        @instance_method: add_overlap

            Adds an overlap edge to the graph.

        @param u: The source read index.

        @param v: The target read index.

        @param d: The direction going from u to v

        @param l: The overhanging suffix length going from u to v

        @raise IndexError: If given invalid read indices.
        """

        if u < 0 or u >= self.n or v < 0 or v >= self.n:
            raise IndexError("Trying to add overlap with out of bounds indices")

        if v in self.g[u] and self.g[u][v][1] <= l:
            return

        self.g[u][v] = (d, l)

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
        dirs = []
        suflens = []
        for u, adj in self.g.items():
            for v, e in adj.items():
                edges.append((u,v))
                dirs.append(e[0])
                suflens.append(e[1])
        G.add_edges(edges)
        G.es['dir'] = dirs
        G.es['len'] = suflens
        G.vs['seq'] = self.seqs
        return G

    #  @classmethod
    #  def generate_overlapped(cls, seqs : tuple, overlap_seeds : list): # -> cls:

        #  n = len(seqs)
        #  g = cls(seqs)

    @classmethod
    def generate_seed_based(cls, seqs : tuple, overlap_seeds : list, k : int): # -> cls:
        """
        @classmethod generate_seed_based:

            Constructs a new OverlapGraph object using k-mer based overlap_seeds
            generated from @FindOverlaps.get_overlap_seeds. Each seeded overlap is
            extended as if the remaining non-seeded parts are assumed to be perfect
            matches. Given multiple seeds between two given reads, The first seed
            parsed is used and the rest are discarded.

        @param seqs: Tuple of reads ordered by read index.

        @param overlap_seeds: List of overlap seeds. See @FindOverlaps.get_overlap_seeds.__doc__
                              for details.

        @param k: The seed length. This is needed for reverse complement overlaps, otherwise
                  we can't find the correct overlaps.

        @return instantiated OverlapGraph object.
        """

        n = len(seqs)
        g = cls(seqs)

        for u, v, upos, vpos, rc in overlap_seeds:

            useq = seqs[u]
            vseq = seqs[v]
            ulen = len(useq)
            vlen = len(vseq)

            if rc: vpos = vlen - vpos - k - 1

            if upos <= vpos and (ulen - upos) <= (vlen - vpos):
                g.add_overlap(v, u, -1, 0)
            elif upos >= vpos and (ulen - upos) >= (vlen - vpos):
                g.add_overlap(u, v, -1, 0)
            else:
                if upos > vpos:
                    suflen = (vlen - vpos) - (ulen - upos)
                    prelen = upos - vpos
                    if not rc:
                        g.add_overlap(u, v, 1, suflen)
                        g.add_overlap(v, u, 2, prelen)
                    else:
                        g.add_overlap(u, v, 0, suflen)
                        g.add_overlap(v, u, 0, prelen)
                else:
                    suflen = vpos - upos
                    prelen = (ulen - upos) - (vlen - vpos)
                    if not rc:
                        g.add_overlap(u, v, 2, suflen)
                        g.add_overlap(v, u, 1, prelen)
                    else:
                        g.add_overlap(u, v, 3, suflen)
                        g.add_overlap(v, u, 3, prelen)
        return g

    @classmethod
    def generate_gold_standard(cls, seqs : tuple, records : list, genome_length : int): # -> cls: (annotations don't work for this use case apparently)
        """
        @classmethod generate_gold_standard:

            Constructs a new OverlapGraph object usingt input reads and
            their gold standard mapping back to original reference.

        @param seqs: Tuple of reads ordered be read index.

        @param records: Read triplets. Look at AssemblyProblem.create_reads.__doc__ function for definition.

        @param genome_length: Reference genome length; used to calculate overlaps
                              that cross the circular barrier.

        @raise ValueError: If len(seqs) != len(records).

        @raise AssertionError: If records are incorrectly sorted (should never happen).

        @return instantiated OverlapGraph object.
        """

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
                    g.add_overlap(u, v, -1, 0)
                else:
                    """
                    If we reach here, then v extends to the right of u by >= 1 nucleotides.
                    If upos == vpos, then u must be contained within v.
                    """
                    if upos == vpos:
                        g.add_overlap(v, u, -1, 0)
                    else:
                        """
                        Reaching here means (upos < vpos) and (upos + ulen < vpos + vlen)
                        """

                        suflen = vpos + vlen - upos - ulen
                        prelen = vpos - upos

                        if not urev:
                            if not vrev:
                                g.add_overlap(u, v, 1, suflen)
                                g.add_overlap(v, u, 2, prelen)
                            else:
                                g.add_overlap(u, v, 0, suflen)
                                g.add_overlap(v, u, 0, prelen)
                        else:
                            if not vrev:
                                g.add_overlap(u, v, 3, suflen)
                                g.add_overlap(v, u, 3, prelen)
                            else:
                                g.add_overlap(u, v, 2, suflen)
                                g.add_overlap(v, u, 1, prelen)

        return g

    def print(self):
        """
        @func print: TODO
        """
        for u in self.g:
            for v in self.g[u]:
                arrow = [' ', '--', ' ']
                d, l = self.g[u][v]
                t, h = ((d>>1)&1, d&1)
                arrow[0] = '>' if t==0 else '<'
                arrow[-1] = '>' if h==1 else '<'
                print("({}) {} ({}) :: {}".format(u, "".join(arrow), v, l))
