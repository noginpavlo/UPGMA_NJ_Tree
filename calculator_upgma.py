import subprocess
import tempfile
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os
import matplotlib.pyplot as plt
from Bio import Phylo
import matplotlib

matplotlib.use('Agg')

class TreeBuilder:
    def __init__(self):
        self.aligned_file = None
        self.dissimilarity_matrix = None
        self.tree = None
        self.n_seq = 2


    def align_sequences(self, list1, list2):
        import os

        # Create a temporary FASTA file to hold the input sequences
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as temp_fasta:
            for seq in list2:
                if isinstance(seq, str):
                    seq = SeqRecord(Seq(seq), id=f"{list1[list2.index(seq)]}")
                temp_fasta.write(f">{seq.id}\n{seq.seq}\n")
            temp_fasta_name = temp_fasta.name  # Save the temp file name

        print(f"Temporary FASTA file created at: {temp_fasta_name}")

        # Run ClustalW to align the sequences
        output_aln = temp_fasta_name.replace('.fasta', '_aligned.aln')  # Define output file
        print(f"Expected output alignment file: {output_aln}")

        clustal_command = [
            '/opt/miniconda3/envs/bioinfo/bin/clustalw',
            f'-INFILE={temp_fasta_name}',
            f'-OUTFILE={output_aln}',
            '-ALIGN'
        ]

        # Run the ClustalW command
        try:
            subprocess.run(clustal_command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running ClustalW: {e}")
            return

        # Read the alignment from the generated .aln file
        try:
            alignment = AlignIO.read(output_aln, "clustal")
        except FileNotFoundError:
            print(f"Alignment file not found: {output_aln}")
            return
        except Exception as e:
            print(f"Error reading alignment file: {e}")
            return

        # Print the aligned sequences
        print("Aligned Sequences:")
        for record in alignment:
            print(f"{record.id}: {record.seq}")

        self.aligned_file = alignment

        # Clean up temporary files
        os.remove(temp_fasta_name)
        os.remove(output_aln)


    def calculate_dissimilarity_matrix(self):
        if self.aligned_file is None:
            print("No alignment data available. Run align_sequences first.")
            return

        print("Loaded alignment:")
        for record in self.aligned_file:
            print(f"{record.id}: {record.seq}")

        calculator = DistanceCalculator("identity")
        distance_matrix = calculator.get_distance(self.aligned_file)

        print("\nDissimilarity Matrix:")
        print(distance_matrix)

        self.dissimilarity_matrix = distance_matrix


    def build_tree(self, method="upgma"):
        if self.dissimilarity_matrix is None:
            print("No distance matrix available. Run calculate_dissimilarity_matrix first.")
            return

        constructor = DistanceTreeConstructor()
        if method.lower() == "upgma":
            tree = constructor.upgma(self.dissimilarity_matrix)
        elif method.lower() == "nj":
            tree = constructor.nj(self.dissimilarity_matrix)
        else:
            raise ValueError(f"Unsupported method: {method}. Use 'upgma' or 'nj'.")

        self.tree = tree

    def visualize_tree(self, output_file=None, output_format="png"):
        if self.tree is None:
            print("No tree data available. Build tree first.")
            return

        """
        Visualizes and optionally exports the phylogenetic tree with custom sample name styles.

        :param output_file: Path to save the tree image (e.g., 'tree.png' or 'tree.jpg')
        :param output_format: Format of the output file ('png', 'jpeg', etc.)
        """
        # Define the directory for saving plots
        save_directory = os.path.join("static", "plots")

        # Ensure the directory exists; create it if it doesn't
        if not os.path.exists(save_directory):
            os.makedirs(save_directory)

        # If no output file is provided, use a default filename
        if output_file is None:
            output_file = os.path.join(save_directory, f"tree.{output_format}")

        # Remove internal node labels (set name to None)
        for clade in self.tree.find_clades():
            if not clade.is_terminal():
                clade.name = None  # Remove the label of internal nodes

        # Create a new Matplotlib figure
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(1, 1, 1)


        Phylo.draw(self.tree, axes=ax, do_show=False,
                   branch_labels=lambda c: f"{c.branch_length:.3f}" if c.branch_length != 0 else None)

        # Now, let's modify the font size and style of the terminal clade labels (only words)
        for text in ax.texts:  # 'ax.texts' contains all the text objects (labels)
            label_text = text.get_text()  # Get the actual text from the label

            if label_text.isalpha():  # Check if the text is a word (no digits)
                text.set_fontsize(24)  # Change font size for words only
                text.set_fontname('Arial')  # Change font style for words only
                text.set_color('blue')  # Change text color for words only

        plt.savefig(output_file, format=output_format, dpi=300)  # High resolution
        print(f"Tree saved to {output_file} in {output_format} format.")

        # Close the figure to release memory and suppress GUI
        plt.close(fig)
