import AssemblyProblem as AP
from OverlapGraph import OverlapGraph

genome = "CGTGAAAAGGTGCTCAGGCG" # len(genome) = 20

#  readinfo = []

#  readinfo.append((0, 10, False))
#  readinfo.append((3, 9, False))
#  readinfo.append((8, 6, False))
#  readinfo.append((10, 7, False))
#  readinfo.append((11, 10, False))

#  seqs, names, records = AP.define_reads(genome, readinfo)

seqs, names, records = AP.create_reads(genome, 5, 9, 2, False)

AP.pretty_layout(seqs, records)

AP.write_fasta("genome.fa", [genome], ["example"])
AP.write_fasta("reads.fa", seqs, names)

overlap_graph = OverlapGraph.generate(seqs, records, len(genome))
overlap_graph.get_igraph().write_gml("overlaps.gml")
