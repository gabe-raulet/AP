import AssemblyProblem as AP
from OverlapGraph import OverlapGraph

genome_length       = 30
read_depth          = 20
mean_read_length    = 12
sd_read_length      = 1
reverse_complements = True

genome = AP.create_random_genome(genome_length)

seqs, names, records = AP.create_reads(genome, read_depth, mean_read_length, sd_read_length, reverse_complements)
AP.pretty_layout(seqs, records, "layout.txt")
AP.pretty_layout(seqs, records)

overlap_graph = OverlapGraph.generate(seqs, records, genome_length)

overlap_igraph = overlap_graph.get_igraph()

AP.write_fasta("genome.fa", [genome], ["chrom1"])
AP.write_fasta("reads.fa", seqs, names)
overlap_igraph.write_gml("overlaps.gml")

G = overlap_graph.get_pruned()
H = G.get_naive_tr(0).get_igraph()
H.delete_vertices(H.vs.select(_degree_eq=0))
H.write_gml("string.gml")

