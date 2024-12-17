from flask import Flask, render_template, request, redirect, url_for
from calculator_upgma import *


app = Flask(__name__)
calculator = TreeBuilder()


# Configure upload folder
UPLOAD_FOLDER = 'static/uploads'  # Directory to save uploaded files
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@app.route("/")
def home():
    return render_template("index2.html")


@app.route("/parameters")
def calculation_parameters():
    return render_template("parameters.html")


@app.route("/input")
def input_seq():
    return render_template("input2.html")


@app.route("/show_result")
def show_result():
    return render_template("result_tree.html")


@app.route('/submit_sequences', methods=['POST'])
def submit_sequences():

    sequence_names = request.form.getlist('sequence_name[]')
    sequence_data = request.form.getlist('sequence_data[]')

    calculator.align_sequences(sequence_names, sequence_data)
    calculator.calculate_dissimilarity_matrix()
    calculator.build_tree()
    calculator.visualize_tree()

    return redirect(url_for('show_result'))


@app.route('/submit_file', methods=['POST'])
def submit_file():

    file = request.files['file-upload']

    # Save the uploaded file
    upload_dir = "static/uploads"  # Directory to save files
    os.makedirs(upload_dir, exist_ok=True)  # Ensure the directory exists
    file_path = os.path.join(upload_dir, file.filename)
    file.save(file_path)

    calculator.align_sequences(fasta_file=file_path)
    calculator.calculate_dissimilarity_matrix()
    calculator.build_tree()
    calculator.visualize_tree()

    return redirect(url_for('show_result'))


@app.route("/elements")
def elements():
    return render_template("elements.html")


@app.route("/inp")
def new_input():
    return render_template("input2.html")


@app.route("/theory")
def theory_section():
    return render_template("theory.html")


@app.route("/docs")
def documentation():
    return render_template("docs.html")


if __name__ == "__main__":
    app.run(debug=True)