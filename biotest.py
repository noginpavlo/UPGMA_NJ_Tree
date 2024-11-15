from Bio import SeqIO


# Step 1: Load sequences
sequences = list(SeqIO.parse("sequences.fasta", "fasta"))