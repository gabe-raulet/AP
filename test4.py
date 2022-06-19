from AssemblyProblem import *
from OverlapGraph import *
from FindOverlaps import *

genome_length       = 1000
read_depth          = 20
mean_read_length    = 50
sd_read_length      = 5
circular            = True
reverse_complements = True

genome = create_random_genome(genome_length)

seqs, names, records = create_reads(genome, read_depth, mean_read_length, sd_read_length, circular, reverse_complements)
#  genome, seqs, names, records = read_problem_instance("genome.fa", "reads.fa")

pretty_layout(seqs, records, "layout.txt")
pretty_layout(seqs, records)

write_fasta("genome.fa", [genome], ["chrom1"])
write_fasta("reads.fa", seqs, names)

#####

k = 9
w = 2
overlap_seeds = get_overlap_seeds(seqs, k, w)

gold_graph = OverlapGraph.generate_gold_standard(seqs, records, genome_length).get_pruned()
seed_graph = OverlapGraph.generate_seed_based(seqs, overlap_seeds, k).get_pruned()

gold_graph.get_igraph().write_gml("gold_graph.gml")
seed_graph.get_igraph().write_gml("seed_graph.gml")

gold_graph.get_naive_tr(0).get_igraph().write_gml("gold_graph_tr.gml")
seed_graph.get_naive_tr(0).get_igraph().write_gml("seed_graph_tr.gml")


