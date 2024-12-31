FROM python:3.12-slim AS build
WORKDIR /cladevo
RUN apt-get update && apt-get install -y \
    clustalw \
    && rm -rf /var/lib/apt/lists/*

# Copy the requirements file into the container
COPY requirements.txt .
# Install the Python dependencies globally
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --no-cache-dir gunicorn

# Runtime stage
FROM python:3.12-slim
WORKDIR /cladevo

# Copy the installed dependencies from the build stage
COPY --from=build /usr/local /usr/local
COPY . .
COPY ssl/certificate.pem /cladevo/
COPY ssl/private_key.pem /cladevo/
EXPOSE 443

# Command to start the application with Gunicorn
CMD ["gunicorn", "-w", "4", "-b", "0.0.0.0:443", "--certfile=certificate.pem", "--keyfile=private_key.pem", "server:app"]

