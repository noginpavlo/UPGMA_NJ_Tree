from flask import Flask, render_template, request, redirect, url_for, jsonify
from calculator_upgma import *

app = Flask(__name__)
calculator = TreeBuilder()


# Configure upload folder
UPLOAD_FOLDER = 'static/uploads'  # Directory to save uploaded files
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

method = "nj"

@app.route("/get-method", methods=['POST'])
def get_method():
    global method
    response = request.form.getlist('method')
    method = response[0]
    print(f"Now method used is {method}")
    return 'Method updated successfully', 200


@app.route("/")
def home():
    return render_template("index2.html")


@app.route("/parameters")
def calculation_parameters():
    return render_template("parameters.html")


@app.route("/input")
def input_seq():
    print(method)
    return render_template("input2.html")


@app.route('/show_result/<unique_id>')
def show_result(unique_id):
    # Log the request URL for debugging purposes
    print(f"Request URL: {request.url}")
    print(f"Received unique_id: {unique_id}")

    # No need to fetch `unique_id` from query parameters; it’s passed in the URL
    if not unique_id:
        print("Error: No unique ID provided")
        return "Error: No unique ID provided", 400

    # Render the template with the unique_id
    return render_template("result_tree.html", unique_id=unique_id)


@app.route('/submit_sequences', methods=['POST'])
def submit_sequences():
    sequence_names = request.form.getlist('sequence_name[]')
    sequence_data = request.form.getlist('sequence_data[]')

    calculator.align_sequences(sequence_names, sequence_data)
    calculator.calculate_dissimilarity_matrix()
    calculator.build_tree(method)

    unique_id = calculator.visualize_tree()

    if unique_id is None:
        return jsonify({'error': 'Tree visualization failed'}), 400

    # Build the redirect URL
    redirect_url = url_for('show_result', unique_id=unique_id)

    # Return the URL as JSON
    return jsonify({'redirect_url': redirect_url})


@app.route('/submit_file', methods=['POST'])
def submit_file():

    file = request.files['file-upload']

    # Save the uploaded file
    upload_dir = "static/uploads"  # Directory to save files
    os.makedirs(upload_dir, exist_ok=True)  # Ensure the directory exists
    file_path = os.path.join(upload_dir, file.filename)
    file.save(file_path)

    # Check file extensions
    if file.filename.lower().endswith('.fasta'):
        # Handle FASTA file
        calculator.align_sequences(fasta_file=file_path)
        calculator.calculate_dissimilarity_matrix()
        calculator.build_tree(method)
        unique_id = calculator.visualize_tree()
    elif file.filename.lower().endswith(('.xls', '.xlsx')):
        # Handle Excel file
        calculator.reformat(file_path)
        calculator.calculate_dissimilarity_matrix()
        calculator.build_tree(method)
        unique_id = calculator.visualize_tree()
    else:
        # Handle unsupported file types (optional)
        return "Unsupported file type", 400

    if unique_id is None:
        return jsonify({'error': 'Tree visualization failed'}), 400

    # Build the redirect URL
    redirect_url = url_for('show_result', unique_id=unique_id)

    # Return the URL as JSON
    return jsonify({'redirect_url': redirect_url})


@app.route("/inp")
def new_input():
    return render_template("input2.html")


@app.route("/theory")
def theory_section():
    return render_template("theory.html")


@app.route("/docs")
def documentation():
    return render_template("docs.html")


@app.route("/test_youtube")
def test_youtube():
    return render_template("test_youtube.html")


if __name__ == "__main__":
    app.run(debug=True)