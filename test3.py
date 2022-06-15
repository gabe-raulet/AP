import AssemblyProblem as AP
from OverlapGraph import OverlapGraph

genome, seqs, names, records = AP.read_problem_instance("genome.fa", "reads.fa")

AP.pretty_layout(seqs, records)
