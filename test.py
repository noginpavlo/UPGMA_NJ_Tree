from calculator_upgma import *

calculator = TreeBuilder()

apple = "AGCTGAAACCCTCT"
banana = "AGCTGCCCTCT"
bee = "AGCTGCCCAAATCT"
watermelon = "AGCTGCCCTCTAAATCT"
cat = "AGCTGCCCTCTAAATCGCTGCTGCTG"
dog = "AGCTGCCCTCTAAACTGCTGCTGCTG"
ship = "AGCTGCCCTCTAAACACACACTGCT"
mouse = "AGCTGCCCTCTAAACTGCTGAACT"
turtle = "AGCTGCCCTCTAAACTGACTGCT"

names = ["apple", "banana", "bee", "watermelon", "cat", "dog", "ship", "mouse", "turtle"]
sequences = ["AGCTGAAACCCTCT", "AGCTGCCCTCT", "AGCTGCCCAAATCT", "AGCTGCCCTCTAAATCT", "AGCTGCCCTCTAAACTGCTGCTGCTG",
             "AGCTGCCCTCTAAACACACACTGCT", "AGCTGCCCTCTAAACTGCTGAACT", "AGCTGCCCTCTAAACTGACTGCT", "ACGTAGCTAGCTAGCACTG"]


calculator.align_sequences(names, sequences)

calculator.calculate_dissimilarity_matrix()

calculator.build_tree()

calculator.visualize_tree()