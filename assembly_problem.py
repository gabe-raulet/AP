import sys
import random
import numpy as np
from igraph import *

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def create_random_genome(size):
    """
        @func create_random_genome:

            Creates a random sequence of nucleotides.

        @param size: the length of the sequence.

        @raise IndexError: if size <= 0.
    """

    if (size <= 0):
        raise IndexError("size <= 0")

    return "".join(random.choice(['A','C','G','T']) for _ in range(size))

def write_fasta(filename, seqs, names):
    """
        @func write_fasta:

            Writes a fasta file.

        @param filename: The name of the fasta to write.

        @param seqs: A list of nucleotide strings.

        @param names: A list of associated fasta record names.

        @raise IndexError: If len(seqs) != lens(names).
    """

    if len(seqs) != len(names):
        raise IndexError("len(seqs) != len(names)")

    with open(filename, "w") as f:
        for i in range(len(seqs)):
            f.write(">{}\n{}\n".format(names[i], seqs[i]))

def circular_slice(s, i, l):
    """
        @func circular_slice:

            Returns a slice from the sequence s, where s is treated
            as a circular sequence. That is, if n = len(s), then semantically
            we have that s[i] == s[i + n] always.

        @param s: The sequence string.

        @param i: The starting index of the slice.

        @param l: The length of the slice.

        @return (cslice, startpos, endpos):

            cslice: The computed circular slice.

            startpos: The starting position of cslice relative to s,
                      i.e. s[startpos] == cslice[0].

            endpos: The modulo reduced ending position (inclusive) of
                    cslice relative to s, i.e. s[endpos] == cslice[-1].
    """

    n = len(s)

    i %= n

    """
        In case i >= len(s), we can just reduce modulo len(s) without
        affecting the result.
    """

    startpos = i

    if i + l <= n:
        """
            If  i + l <= len(s), then s[i : i + l] is a valid python slice.
        """

        endpos = i + l - 1

        cslice = s[i : endpos + 1]

    else:
        """
            If i + l > len(s), then we take the python slice s[i : len(s)] as prefix
            and s[0 : i + l - len(s)] as the wrap-around suffix.
        """

        endpos = i + l - n - 1

        cslice = s[i : ] + s[ : endpos + 1]

    return (cslice, startpos, endpos)

def create_reads(genome, read_depth, mean_read_length, sd_read_length, reverse_complements=True):

    """
        @func create_reads:

            Simulates random perfect reads from the sequence @genome, which is treated as a circular
            sequence. Read lengths are drawn from a normal distribution (not necessarily realistic but
            fine for now).

        @param genome: The genome sequence string.

        @param read_depth: The average read depth.

        @param mean_read_length: The average read length.

        @param sd_read_length: The standard deviation of the read length normal distribution.

        @param reverse_complements: Boolean for whether reverse complements should be simulated
                                    with 1/2 probability.

        @return (seqs, names, records):

            seqs: The raw simulated read sequence strings, ordered by index (immutable tuple).

            names: The generated read "names" (think more like descriptions),
                   ordered by index (immutable tuple).

            records: Triplets recording basic info for each read sequence. For example, for read
                     i, records[i] = (i, startpos, readrev), where startpos is the starting
                     position of read i relative to the genome, and readrev is a boolean
                     representing whether the read is reverse complemented relative to the genome.
                     Returned as a list because it is intended for this list to be sorted by
                     start position.
    """

    genome_length = len(genome)
    num_reads = int((genome_length * read_depth) / mean_read_length)

    seqs = []
    names = []
    records = []

    for i in range(num_reads):

        readpos = random.randint(0,genome_length-1)

        while True:
            readlen = int(np.random.normal(mean_read_length, sd_read_length))
            if readlen > 0:
                break

        readrev = False if not reverse_complements else bool(random.randint(0,1))

        readseq, startpos, endpos = circular_slice(genome, readpos, readlen)

        if startpos < endpos:
            coords = "[{}..{}]".format(startpos, endpos)
        else:
            coords = "[{}..) ++ [..{}]".format(startpos, endpos)

        readname = "R{} | coords :: {} | length :: {} | rev :: {}".format(i, coords, readlen, readrev)

        if readrev:
            readseq = readseq.translate(comp_tab)[::-1]

        seqs.append(readseq)
        names.append(readname)
        records.append((i, startpos, readrev))

    return tuple(seqs), tuple(names), list(records)


def pretty_print_layout(seqs, records):
    """
        @func pretty_print_layout:

            Print out sequences according to their gold standard
            layout.

        @param seqs: The sequences.

        @param records: The sequence record triplets (defined at @create_reads).

        @raise ValueError: len(seqs) != len(records).
    """

    n = len(seqs)

    if n != len(records):
        raise ValueError("len(seqs) != len(records).")

    records = sorted(records, key=lambda t: t[1])

    for i in range(n):
        u, upos, _ = records[i]
        useq = seqs[u]
        print("{:>8}: {}{}".format(u, ' ' * upos, useq))
