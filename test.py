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


calculator.align_sequences(apple, banana, bee, watermelon, cat, dog, ship, turtle, mouse)

calculator.calculate_dissimilarity_matrix()

calculator.build_tree()

calculator.visualize_tree()