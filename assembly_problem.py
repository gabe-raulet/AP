import sys
import random
import numpy as np
from igraph import *

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def create_random_genome(size):
    """
        Creates a random sequence of nucleotides.

        @param size: the length of the sequence

        @raise IndexError: if size <= 0
    """

    if (size <= 0):
        raise IndexError("size <= 0")

    return "".join(random.choice(['A','C','G','T']) for _ in range(size))

def write_fasta(filename, seqs, names):
    """
        Writes a fasta file.

        @param filename: the name of the fasta to write

        @param seqs: a list of nucleotide strings

        @param names: a list of associated fasta record names

        @raise IndexError: if len(seqs) != lens(names)
    """

    if len(seqs) != len(names):
        raise IndexError("len(seqs) != len(names)")

    with open(filename, "w") as f:
        for i in range(len(seqs)):
            f.write(">{}\n{}\n".format(names[i], seqs[i]))

def circular_slice(s, i, l):
    n = len(s)
    if i <= n - l: return (s[i:i+l], i, i+l-1)
    else: return (s[i:] + s[:i+l-G], i, i+l-G-1)

def create_reads(genome, read_depth, mean_read_length, sd_read_length, reverse_complements=True):

    genome_length = len(genome)
    num_reads = int((genome_length * read_depth) / mean_read_length)

    seqs = []

    for i in range(num_reads):

        readpos = random.randint(0,genome_length-1):

            while True:
                readlen = int(np.random.normal(mean_read_length, sd_read_length))
                if readlen > 0:
                    break

        readrev = False if not reverse_complements else bool(random.randint(0,1))

        readseq, startpos, endpos = circular_slice(genome, startpos, endpos)

        if startpos < endpos:
            coords = "[{}..{}]".format(startpos, endpos)
        else:
            coords = "[{}..) ++ [..{}]".format(startpos, endpos)

        readname = "R{} | coords :: {} | length :: {} | rev :: {}".format(i+1, coords, readlen, readrev)

    return seqs
