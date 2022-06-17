from AssemblyProblem import reverse_complement, read_fasta, comp_tab

nt4char = ['A', 'C', 'G', 'T']
nt4map = {nt4char[i] : i for i in range(4)}

def kmers(s : str, k : int):
    """
    @generator kmers:

        Returns a generator yielding every kmer in s.

    @param s: The input string.

    @param k: The k-mer length.

    @yields (i, s[i..i+k)) for each position i in s where
            there is a k-mer.
    """
    l = len(s)
    for i in range(l-k+1):
        yield i, s[i:i+k]

def minimizers(s : str, k : int, w : int):
    """
    @generator minimizers:

        Returns a generator yielding every minimizer in s of
        length k within a window of w.

    @param s: The input string.

    @param k: The minimizer length.

    @param w: The selection window length.

    @yields (i, s[i..i+k) for each position i in s where there
            is a minimizer.
    """

    l = len(s)
    kmers = set()
    for i in range(l-k-w+1):
        pos, minimizer = min([(i+j, s[i+j:i+j+k]) for j in range(w)], key=lambda x: x[1])
        if not minimizer in kmers:
            kmers.add(minimizer)
            yield pos, minimizer

def kmercode(s : str) -> tuple:
    """
    @func kmercode:

        Calculates the integer value code of a k-mer. Uses the
        canonical form a k-mer, where the canonical form of any
        DNA sequence is defined as the smaller of itself or its
        reverse complement.

    @param s: The input k-mer string.

    @return (code, revflag):

        code: The integer value code. This code is unique for every
              every k-mer.

        revflag: False if the canonical(s) == s,
                 True if canonical(s) = reverse_complement(s)
    """
    k = len(s)
    r = s.translate(comp_tab)
    forward = sum(nt4map[s[i]] * (4**(k-i-1)) for i in range(k))
    reverse = sum(nt4map[r[k-i-1]] * (4**(k-i-1)) for i in range(k))
    return min(forward, reverse), (reverse < forward)

def codekmer(code : int, k : int) -> str:
    """
    @func codekmer:

        Calculates the k-mer string of a given k-mer integer value code.
        Note that if the code of a non-canonical sequence is passed, then
        codekmer will return the non-canonical sequence. In order to
        avoid this, you should pass canonical codes. Note that the inverse
        function @kmercode (defined above) always returns the canonical code
        of a sequence.

    @param code: The integer value code.

    @param k: The k-mer length.

    @return s: The k-mer.
    """
    s = [0 for _ in range(k)]
    for i in range(k):
        s[k-i-1] = code % 4
        code //= 4
    return "".join(nt4char[c] for c in s)

def get_kmer_array(seqs : list, k : int, w : int) -> list:
    """
    @func get_kmer_array:

        Construct an array of k-mers (really minimizer k-mers) from an
        input list of DNA strings.

    @param seqs: The list of DNA strings to count k-mers from.

    @param k: The k-mer length.

    @param w: The minimizer window length.

    @return kmer_array:

        Each entry in the kmer_array is a tuple of the form

            (code, readid, readpos, revflag)

        where:

            code: The canonical k-mer code.

            readid: The readid from which the k-mer was obtained.

            readpos: The offset within the read at which the k-mer
                     starts.

            revflag: False if the actual k-mer parsed is the canonical k-mer;
                     True if the actual k-mer parsed is the non-canonical k-mer.

        For example, if revflag == False then seqs[readid][readpos:readpos+k] is
        the k-mer associated with code, and if revflag = True, then the
        reverse complement of the previous expression is the k-mer associated
        with code.
    """
    n = len(seqs)
    kmer_array = []

    for read_index, readseq in enumerate(seqs):
        l = len(readseq)
        if l < k: continue
        for pos, kmer in minimizers(readseq, k, w):
            code, rev = kmercode(kmer)
            kmer_array.append((code, read_index, pos, rev))

    return kmer_array

def compress_kmer_array(kmer_array : list) -> dict:
    """
    @func compress_kmer_array:

        Takes a k-mer array (list) and compresses it into
        a dictionary.

    @param kmer_array: A list of (code, readid, readpos, revflag) tuples, presumably
                       obtained from @get_kmer_array.

    @return kmer_dict:

        Each key in kmer_dict is a k-mer code, and each value is a list of
        (readid, readpos, revflag) tuples, where the compression comes from
        taking all entries with the same k-mer code in kmer_array and making
        them accessible by calling kmer_dict[code].
    """
    A = sorted(kmer_array)
    kmer_dict = {}
    i = 0
    while i < len(A)-1:
        code = A[i][0]
        adj = []
        for j in range(i, len(A)):
            if A[j][0] != code:
                break
            adj.append(A[j][1:])
        kmer_dict[code] = adj
        i = j
    return kmer_dict

def get_overlap_seeds(seqs : list, k : int, w : int) -> list:
    """
    @func get_overlap_seeds:

        Takes a list of DNA string sequences and uses minimizers as
        seeds in order to compute a list of every location where
        two reads share the same k-mer.

    @param seqs: The DNA string sequences.

    @param k: The minimizer length.

    @param w: The minimizer window length.

    @return overlap_seeds:

        A list of (u, v, upos, vpos, urev, vrev) tuples where:

            u: The source read index.

            v: The target read index.

            upos: The starting position within u of the shared k-mer seed.

            vpos: The starting position within v of the shared k-mer seed.

            urev: Whether the k-mer substring in u is non-canonical.

            vrev: Whether the k-mer substring in v is non-canonical.
    """

    kmer_dict = compress_kmer_array(get_kmer_array(seqs, k, w))

    overlap_seeds = []

    for code, triples in kmer_dict.items():
        l = len(triples)
        for j in range(1, l):
            for i in range(j):
                u, upos, urev = triples[i]
                v, vpos, vrev = triples[j]
                overlap_seeds.append((u, v, upos, vpos, urev, vrev))

    return overlap_seeds

