from AssemblyProblem import *
from OverlapGraph import *
from FindOverlaps import *

genome_length       = 15
read_depth          = 5
mean_read_length    = 6
sd_read_length      = 1
circular            = False
reverse_complements = True

genome = create_random_genome(genome_length)

#  seqs, names, records = create_reads(genome, read_depth, mean_read_length, sd_read_length, circular, reverse_complements)
genome, seqs, names, records = read_problem_instance("genome.fa", "reads.fa")

pretty_layout(seqs, records, "layout.txt")
pretty_layout(seqs, records)

write_fasta("genome.fa", [genome], ["chrom1"])
write_fasta("reads.fa", seqs, names)

#####

k = 4
w = 0
overlap_seeds = get_overlap_seeds(seqs, k, w)
