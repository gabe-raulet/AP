import sys
import random
import numpy as np
from igraph import *

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def reverse_complement(s : str) -> str:
    return s.translate(comp_tab)[::-1]

def create_random_genome(size : int) -> str:
    """
    @func create_random_genome:

        Creates a random sequence of nucleotides.

    @param size: the length of the sequence.

    @raise IndexError: if size <= 0.
    """

    if (size <= 0):
        raise IndexError("size <= 0")

    return "".join(random.choice(['A','C','G','T']) for _ in range(size))

def write_fasta(filename : str, seqs : list, names : list) -> None:
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

def circular_slice(s : str, i : int, l : int) -> tuple:
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


def create_reads(genome : str, read_depth : int, mean_read_length : int, sd_read_length : float, circular : bool = True, reverse_complements : bool = True) -> tuple:

    """
    @func create_reads:

        Simulates random perfect reads from the sequence @genome, which is treated as a circular
        sequence. Read lengths are drawn from a normal distribution (not necessarily realistic but
        fine for now).

    @param genome: The genome sequence string.

    @param read_depth: The average read depth.

    @param mean_read_length: The average read length.

    @param sd_read_length: The standard deviation of the read length normal distribution.

    @param circular: Boolean for whether genome should be treated as circular or not

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

        if circular:
            readpos = random.randint(0,genome_length-1)
            while True:
                readlen = int(np.random.normal(mean_read_length, sd_read_length))
                if readlen > 0:
                    break
            readseq, startpos, endpos = circular_slice(genome, readpos, readlen)

        else:
            readpos = random.randint(0,genome_length-mean_read_length-1)
            while True:
                readlen = int(np.random.normal(mean_read_length, sd_read_length))
                if readlen > 0 and readpos + readlen <= genome_length:
                    break

            startpos = readpos
            endpos = readpos + readlen - 1
            readseq = genome[readpos : endpos + 1]

        readrev = False if not reverse_complements else bool(random.randint(0,1))

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

def define_reads(genome : str, readinfo : list) -> tuple:

    """
    @func define_reads:

        Creates perfect reads from the sequence @genome (treated as circular sequence) as specified
        by the user input.

    @param genome: The genome sequence string.

    @param readinfo: Triplets (readpos, readlen, readrev) specifying the starting position of
                     each read relative to the genome, the length of the read, and whether
                     the read has been reverse complemented relative to the genome.

    @return (seqs, names, records): Defined in create_reads.__doc__.

    """

    genome_length = len(genome)
    num_reads = len(readinfo)

    seqs = []
    names = []
    records = []

    for i in range(num_reads):

        readpos, readlen, readrev = readinfo[i]

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

def pretty_layout(seqs : tuple, records : list, filename = None) -> None:
    """
    @func pretty_layout:

        Prettily layout sequences according to their gold standard
        layout.

    @param seqs: The sequences.

    @param records: The sequence record triplets (defined in create_reads.__doc__).

    @param filename: The (optional) output file path. If filename is a str object,
                     then a file will be created with the output redirected there.
                     Otherwise, output will be written to stdout.

    @raise ValueError: len(seqs) != len(records).
    """

    n = len(seqs)

    if n != len(records):
        raise ValueError("len(seqs) != len(records).")

    records = sorted(records, key=lambda t: t[1])

    if isinstance(filename, str):
        with open(filename, "w") as f:
            for i in range(n):
                u, upos, urev = records[i]
                useq = seqs[u]
                uformat = ' ' * (upos) + '<' + useq.translate(comp_tab)[::-1] if urev else ' ' * (upos+1) + useq + '>'
                f.write("{:>4}: {}\n".format(u, uformat))
    else:
        for i in range(n):
            u, upos, urev = records[i]
            useq = seqs[u]
            uformat = ' ' * (upos) + '<' + useq.translate(comp_tab)[::-1] if urev else ' ' * (upos+1) + useq + '>'
            sys.stdout.write("{:>4}: {}\n".format(u, uformat))
        sys.stdout.flush()

def read_fasta(fasta : str) -> tuple:
    """
    @func read_fasta:

        Reads the sequences and names from a FASTA file.

    @param fasta: FASTA filename.

    @raise AssertionError: If len(seqs) != len(names).

    @return (seqs, names):

        seqs: The raw read sequences.

        names: The corresponding read sequence names.
    """
    seqs = []
    names = []
    with open(fasta, "r") as f:
        seqlines = []
        name = ""
        for line in f.readlines():
            if line.startswith('>'):
                if len(seqlines) > 0:
                    seqs.append("".join(seqlines))
                    names.append(name)
                    seqlines = []
                name = line.lstrip('>').rstrip()
            else:
                seqlines.append(line.rstrip())
        if len(seqlines) > 0:
            seqs.append("".join(seqlines))
            names.append(name)
    assert len(seqs) == len(names)
    return seqs, names

def read_problem_instance(genome_filename : str, reads_filename : str) -> tuple:
    """
    @func read_problem_instance:

        Reads a problem instance from disk. Assumes problem instance
        was previously created by this module and written to disk by it.

    @param genome_filename: The FASTA file corresponding to the reference genome.

    @param reads_filename: The FASTA file corresponding to the reads.

    @raise ValueError: If the genome sequence and the input reads are incompatible.
                       For each read sequence read from disk, its coordinate information
                       is cross-referenced with the genome to make sure the two files
                       are consistent and refer to the same problem instance.

                       Additionally, we require that the genome has exactly one sequence/chromosome.

    @return (genome, seqs, names, records):

        genome: The genome sequence string.

        seqs: The ordered tuple of read sequence strings.

        names: The ordered tuple of read sequence names (not vital when
               reading a problem from disk but here for consistency with
               other functions).

        records: The sequence record triplets (defined in create_reads.__doc__).
    """

    _genome, _genome_name = read_fasta(genome_filename)

    if len(_genome) != 1:
        raise ValueError("Genome FASTA has more than one sequence")

    genome = _genome.pop()
    genome_name = _genome_name.pop()

    seqs, names = read_fasta(reads_filename)

    records = []

    for i in range(len(names)):
        readseq = seqs[i]
        readlen = len(readseq)
        name = names[i]
        readpos = int(name.split('|')[1].split('[')[1].split('..')[0])
        readrev = bool(name.split("::")[-1].lstrip().rstrip() == "True")

        readseq_check, startpos_check, endpos_check = circular_slice(genome, readpos, readlen)

        if readrev:
            readseq_check = readseq_check.translate(comp_tab)[::-1]

        if readseq_check != readseq or startpos_check != readpos:
            raise ValueError("Genome fasta and reads fasta are incompatible")

        records.append((i, readpos, readrev))

    return str(genome), tuple(seqs), tuple(names), list(records)
