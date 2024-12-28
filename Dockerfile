FROM python:3.12-slim

# Set the working directory inside the container
WORKDIR /cladevo

# Copy the requirements file into the container
COPY requirements.txt /cladevo/

# Install dependencies from requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install gunicorn explicitly
RUN pip install gunicorn

# Copy the rest of the application code into the container
COPY . /cladevo/

# Expose the port the app will run on
EXPOSE 5000

# Run the Flask app using Gunicorn
CMD ["gunicorn", "-w", "4", "-b", "0.0.0.0:5000", "server:app"]

