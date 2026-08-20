FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies + gosu
RUN apt-get update && \
    apt-get install -y build-essential gfortran sudo git wget gosu && \
    rm -rf /var/lib/apt/lists/*

# Create user 'lmto'
RUN useradd -m -s /bin/bash lmto && \
    echo "lmto:lmto" | chpasswd && \
    adduser lmto sudo

# Create calculations directory owned by lmto
RUN mkdir -p /home/lmto/lmto_calculations && \
    chown -R lmto:lmto /home/lmto/lmto_calculations

USER lmto
WORKDIR /home/lmto

# Build LMTO binary
RUN mkdir -p bin
COPY --chown=lmto:lmto source.tar.gz bin/source.tar.gz
RUN cd bin && \
    tar -xzf source.tar.gz --strip-components=1 && \
    rm source.tar.gz && \
    make all
RUN echo 'export PATH="$HOME/bin:$PATH"' >> .bashrc

# Install Miniforge
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p /home/lmto/miniconda && \
    rm miniforge.sh

ENV PATH="/home/lmto/miniconda/bin:$PATH"
RUN conda init bash && \
    conda create -n lmto_env python=3.12 -y && \
    echo "conda activate lmto_env" >> .bashrc

# Install htlmto
RUN git clone https://github.com/balaranjan/High-throughput-LMTO.git && \
    cd High-throughput-LMTO && \
    /home/lmto/miniconda/envs/lmto_env/bin/pip install .

# Entrypoint runs as root, fixes permissions, then drops to lmto
USER root
RUN printf '#!/bin/bash\nchown -R lmto:lmto /home/lmto/lmto_calculations\nexec gosu lmto /bin/bash --login\n' \
    > /entrypoint.sh && chmod +x /entrypoint.sh

WORKDIR /home/lmto/lmto_calculations
ENTRYPOINT ["/entrypoint.sh"]
