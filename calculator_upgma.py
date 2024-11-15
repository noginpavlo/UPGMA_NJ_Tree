import subprocess
import tempfile
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt


class TreeBuilder:
    def __init__(self):
        self.aligned_file = None
        self.dissimilarity_matrix = None
        self.tree = None


    def align_sequences(self, *args):
        import os

        # Step 1: Create a temporary FASTA file to hold the input sequences
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as temp_fasta:
            for idx, seq in enumerate(args):
                if isinstance(seq, str):
                    seq = SeqRecord(Seq(seq), id=f"seq{idx}")
                temp_fasta.write(f">{seq.id}\n{seq.seq}\n")
            temp_fasta_name = temp_fasta.name  # Save the temp file name

        print(f"Temporary FASTA file created at: {temp_fasta_name}")

        # Step 2: Run ClustalW to align the sequences
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

        # Step 3: Read the alignment from the generated .aln file
        try:
            alignment = AlignIO.read(output_aln, "clustal")
        except FileNotFoundError:
            print(f"Alignment file not found: {output_aln}")
            return
        except Exception as e:
            print(f"Error reading alignment file: {e}")
            return

        # Step 4: Print the aligned sequences
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

        return distance_matrix

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

        # Visualize the tree
        print("\nPhylogenetic Tree:")
        Phylo.draw_ascii(tree)

        self.tree = tree

        return tree

    def visualize_tree(self, output_file=None, output_format="png"):
        if self.tree is None:
            print("No tree data available. Build tree first.")
            return

        """
        Visualizes and optionally exports the phylogenetic tree.

        :param output_file: Path to save the tree image (e.g., 'tree.png' or 'tree.jpg')
        :param output_format: Format of the output file ('png', 'jpeg', etc.)
        """
        # Visualize the tree
        fig = plt.figure(figsize=(10, 8))  # Adjust size as needed
        ax = fig.add_subplot(1, 1, 1)

        Phylo.draw(self.tree, axes=ax)

        # Export to file if output_file is specified
        if output_file:
            plt.savefig(output_file, format=output_format, dpi=300)  # High resolution
            print(f"Tree saved to {output_file} in {output_format} format.")

        # Show the plot (optional)
        plt.show()