import sys
import AssemblyProblem as AP
from OverlapGraph import OverlapGraph

reverse_complements = False

def main(reads_filename : str, genome_filename : str, gml_prefix : str, genome_length : int, read_depth : int, mean_read_length : int, sd_read_length : float) -> None:

    genome = AP.create_random_genome(genome_length)

    seqs, names, records = AP.create_reads(genome, read_depth, mean_read_length, sd_read_length, reverse_complements)

    AP.pretty_layout(seqs, records)

    overlap_graph = OverlapGraph.generate_gold_standard(seqs, records, genome_length)

    overlap_igraph = overlap_graph.get_igraph()

    AP.write_fasta(genome_filename, [genome], ["chrom1"])

    AP.write_fasta(reads_filename, seqs, names)

    overlap_igraph.write_gml("{}.overlaps.gml".format(gml_prefix))

if __name__ == "__main__":
    if len(sys.argv) != 8:
        exe = sys.argv[0]
        sys.stderr.write("usage: {} <reads.fa> <ref.fa> <ava> <G> <D> <L> <mu>\n\n".format(exe))
        sys.stderr.write("    <reads.fa> STR :: IN/OUT :: fasta file of reads\n")
        sys.stderr.write("    <ref.fa>   STR :: IN/OUT :: fasta file of random genome\n")
        sys.stderr.write("    <ava>      STR :: IN/OUT :: prefix for output gml files\n")
        sys.stderr.write("    <G>        INT :: IN     :: length of random genome\n")
        sys.stderr.write("    <D>        INT :: IN     :: average read depth\n")
        sys.stderr.write("    <L>        INT :: IN     :: average read length\n")
        sys.stderr.write("    <mu>     FLOAT :: IN     :: standard deviation of read length\n")
        sys.exit(-1)


    reads_filename = sys.argv[1]
    genome_filename = sys.argv[2]
    gml_prefix = sys.argv[3]
    genome_length = int(sys.argv[4])
    read_depth = int(sys.argv[5])
    mean_read_length = int(sys.argv[6])
    sd_read_length = float(sys.argv[7])

    main(reads_filename, genome_filename, gml_prefix, genome_length, read_depth, mean_read_length, sd_read_length)

    sys.exit(0)

