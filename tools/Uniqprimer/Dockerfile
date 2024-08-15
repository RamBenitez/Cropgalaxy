# Use the official Python image from the Docker Hub
FROM python:3.9-slim

# Set the working directory inside the container
WORKDIR /app

# Copy all files from the current directory to /app in the container
COPY . /app

# Make the shell script executable
RUN chmod +x /app/uniqprimer.sh

# Install the required Python packages
RUN pip install --no-cache-dir -r /app/uniqprimer-0.5.0/requirements.txt

# Set the entry point to the shell script
ENTRYPOINT ["/app/uniqprimer.sh"]
