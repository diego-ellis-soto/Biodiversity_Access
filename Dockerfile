FROM quay.io/jupyter/minimal-notebook:ubuntu-24.04

USER root

# Install R
RUN curl -s https://raw.githubusercontent.com/boettiger-lab/repo2docker-r/refs/heads/main/install_r.sh | bash

WORKDIR /code

# Install R packages first.
# This layer only rebuilds when install.r changes.
COPY install.r .
RUN Rscript install.r

# Copy application and data after packages are installed.
COPY . .

# Hugging Face Spaces
EXPOSE 7860

CMD ["R", "--quiet", "-e", "shiny::runApp('app.R', host='0.0.0.0', port=7860)"]
