UniqPrimer
Overview
UniqPrimer is a tool designed for the specific design of PCR primers targeting bacterial sequences. The tool is intended to process "include" FASTA files of diagnostic target genomes (either draft or complete) and "exclude" files of non-target genomes. The primary goal of UniqPrimer is to generate a list of PCR primers that will amplify the target genomes while avoiding non-target genomes.
Features
Design Specific Primers: Generates primers for target genomes while excluding non-target genomes.
FASTA Input: Accepts FASTA files for both target and non-target genomes.
PCR Primer Output: Provides a list of PCR primers suitable for the specified targets.
Development and Encapsulation
UniqPrimer was developed by the Jan Leach Lab at Colorado State University. It has been encapsulated into the Galaxy platform by IRRI (International Rice Research Institute) and the South Green bioinformatics platform.
License
UniqPrimer is licensed under the GNU General Public License v3.
Dockerized UniqPrimer
UniqPrimer is available as a Docker container for easy deployment and execution. The Docker container ensures that UniqPrimer runs in a consistent environment regardless of the host system's configuration.
Building the Docker Image
To build the Docker image for UniqPrimer, use the following command:

docker build -t uniqprimer .
Running UniqPrimer in Docker
To run UniqPrimer with Docker, use the following command:
docker run --rm uniqprimer ./uniqprimer.sh --help

This command will execute UniqPrimer and display the help information.
Docker Compose Setup
If you prefer using Docker Compose, you can use the provided docker-compose.yml file. This file defines how the UniqPrimer container should be built and run.
To start UniqPrimer with Docker Compose, use:
docker-compose build
docker-compose up
This will build and start the UniqPrimer container. If you need to pass specific arguments or commands, you can adjust the docker-compose.yml file accordingly.
Example
To run UniqPrimer with your data files, you can use a command similar to the following
docker run --rm -v /path/to/data:/data uniqprimer ./uniqprimer.sh --include /data/include.fasta --exclude /data/exclude.fasta --output /data/output.txt

This command mounts a local directory (/path/to/data) into the container and specifies the input and output files.


Links
Installing Docker Guide: https://docs.google.com/document/d/1hm500lgpB08p_LezRNM9X12imP3ssYdca87Q4K02hls/edit?usp=sharing
Github Link: https://github.com/SouthGreenPlatform/Uniqprimer.git
