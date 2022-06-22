import AssemblyProblem as AP
from OverlapGraph import OverlapGraph
from FindOverlaps import *

genome_length       = 100000
read_depth          = 30
mean_read_length    = 14050
sd_read_length      = 1000
reverse_complements = True
circular            = True

genome = AP.create_random_genome(genome_length)
seqs, names, records = AP.create_reads(genome, read_depth, mean_read_length, sd_read_length, circular, reverse_complements)
#  AP.pretty_layout(seqs, records, "layout.txt")
#  AP.pretty_layout(seqs, records)

overlap_graph = OverlapGraph.generate_gold_standard(seqs, records, genome_length)
AP.write_fasta("AP_100k/genome.fa", [genome], ["chrom1"])
AP.write_fasta("AP_100k/reads.fa", seqs, names)
overlap_graph.get_igraph().write_gml("AP_100k/gold_overlap.gml")
G = overlap_graph.get_pruned()
G.get_naive_tr(0).get_igraph().write_gml(("AP_100k/gold_string.gml"))

overlap_seeds = get_overlap_seeds(seqs, 31, 19)
dirty_graph = OverlapGraph.generate_seed_based(seqs, overlap_seeds, 31)
dirty_graph.get_igraph().write_gml("AP_100k/dirty_overlap.gml")
H = dirty_graph.get_pruned()
H.get_naive_tr(0).get_igraph().write_gml("AP_100k/dirty_string.gml")
