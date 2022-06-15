import AssemblyProblem as AP
from OverlapGraph import OverlapGraph

genome_length       = 20
read_depth          = 10
mean_read_length    = 10
sd_read_length      = 4
reverse_complements = False

genome = AP.create_random_genome(genome_length)

seqs, names, records = AP.create_reads(genome, read_depth, mean_read_length, sd_read_length, reverse_complements)
AP.pretty_layout(seqs, records, "layout.txt")

overlap_graph = OverlapGraph.generate(seqs, records, genome_length)

overlap_igraph = overlap_graph.get_igraph()

AP.write_fasta("genome.fa", [genome], ["chrom1"])
AP.write_fasta("reads.fa", seqs, names)
overlap_igraph.write_gml("overlaps.gml")

