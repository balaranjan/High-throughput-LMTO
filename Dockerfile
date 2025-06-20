FROM ubuntu:20.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && \
    apt-get install -y build-essential gfortran sudo git wget && \
    rm -rf /var/lib/apt/lists/*

# Create user 'lmto' with password 'lmto'
RUN useradd -m -s /bin/bash lmto && \
    echo "lmto:lmto" | chpasswd && \
    adduser lmto sudo

# Create data exchange directory (accessible via bind mount)
RUN mkdir -p /home/lmto/lmto_calculations && \
    chown -R lmto:lmto /home/lmto/lmto_calculations

# Switch to user context
USER lmto
WORKDIR /home/lmto

# Create bin directory for source code
RUN mkdir -p bin

# Copy and extract source code (replace USER_DEFINED_PATH with actual path)
COPY --chown=lmto:lmto source.tar.gz bin/source.tar.gz
RUN cd bin && \
    tar -xzf source.tar.gz --strip-components=1 && \
    rm source.tar.gz && \
    make all

# Add bin to PATH in .bashrc
RUN echo 'export PATH="$HOME/bin:$PATH"' >> .bashrc

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /home/lmto/miniconda && \
    rm miniconda.sh

# Initialize conda and setup environment
ENV PATH="/home/lmto/miniconda/bin:$PATH"
RUN conda init bash && \
    conda create -n lmto_env python=3.12 -y && \
    echo "conda activate lmto_env" >> .bashrc

# Clone repository and install package
RUN git clone https://github.com/balaranjan/High-throughput-LMTO.git && \
    cd High-throughput-LMTO && \
    /home/lmto/miniconda/envs/lmto_env/bin/pip install .

# Set default command
CMD ["/bin/bash"]
WORKDIR /home/lmto/lmto_calculations
