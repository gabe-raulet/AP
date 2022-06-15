import AssemblyProblem as AP
from OverlapGraph import OverlapGraph

genome = "CGTGAAAACGTGCTCTGGAG" # len(genome) = 20

readinfo = []

readinfo.append((0, 13, False))
readinfo.append((5, 13, False))
readinfo.append((10, 13, False))
readinfo.append((15, 13, False))

seqs, names, records = AP.define_reads(genome, readinfo)

AP.pretty_layout(seqs, records)

AP.write_fasta("genome.fa", [genome], ["example"])
AP.write_fasta("reads.fa", seqs, names)

overlap_graph = OverlapGraph.generate(seqs, records, len(genome))
overlap_graph.get_igraph().write_gml("overlaps.gml")
