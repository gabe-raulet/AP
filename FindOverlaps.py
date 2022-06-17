from AssemblyProblem import reverse_complement, read_fasta, comp_tab

nt4char = ['A', 'C', 'G', 'T']
nt4map = {nt4char[i] : i for i in range(4)}

def canonical(s : str) -> str:
    return min(s, reverse_complement(s))

def minimizers(s : str, k : int, w : int):
    l = len(s)
    kmers = set()
    for i in range(l-k-w+1):
        pos, minimizer = min([(i+j, s[i+j:i+j+k]) for j in range(w)], key=lambda x: x[1])
        if not minimizer in kmers:
            kmers.add(minimizer)
            yield pos, minimizer

def kmers(s : str, k : int):
    l = len(s)
    for i in range(l-k+1):
        yield i, s[i:i+k]

def kmercode(s : str) -> tuple:
    k = len(s)
    r = s.translate(comp_tab)
    forward = sum(nt4map[s[i]] * (4**(k-i-1)) for i in range(k))
    reverse = sum(nt4map[r[k-i-1]] * (4**(k-i-1)) for i in range(k))
    return min(forward, reverse), (reverse < forward)

def codekmer(code : int, k : int) -> str:
    s = [0 for _ in range(k)]
    for i in range(k):
        s[k-i-1] = code % 4
        code //= 4
    return "".join(nt4char[c] for c in s)

def get_kmer_array(seqs : list, k : int, w : int) -> list:
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


#  def get_overlap_seeds(seqs : list, k : int, w : int) -> list:

    #  kmer_dict = compress_kmer_array(get_kmer_array(seqs, k, w))

