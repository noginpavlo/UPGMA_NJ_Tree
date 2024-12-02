from flask import Flask, render_template, request, redirect, url_for
from calculator_upgma import *


app = Flask(__name__)
calculator = TreeBuilder()


@app.route("/")
def home():
    return render_template("index2.html")


@app.route("/parameters")
def calculation_parameters():
    return render_template("parameters.html")


@app.route("/input")
def input_seq():
    return render_template("input.html")


@app.route("/result")
def show_result():
    print("Was called")
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


@app.route("/elements")
def elements():
    return render_template("elements.html")


@app.route("/inp")
def new_input():
    return render_template("input2.html")


@app.route("/button")
def button():
    return render_template("button.html")


if __name__ == "__main__":
    app.run(debug=True)