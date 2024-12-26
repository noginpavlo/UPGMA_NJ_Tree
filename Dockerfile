# Use an official Python runtime as the base image
FROM python:3.12-slim

# Set the working directory inside the container
WORKDIR /cladevo

# Copy the requirements file into the container
COPY requirements.txt /cladevo/

# Install the dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application code into the container
COPY . /cladevo/

# Set the environment variable for Flask to point to server.py
ENV FLASK_APP=server.py

# Expose the port the app will run on
EXPOSE 5000

# Run the Flask app
CMD ["flask", "run", "--host=0.0.0.0"]

