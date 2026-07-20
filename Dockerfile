# syntax=docker/dockerfile:1

# ============================================================
# microscape-nf container
# Pre-installs both Python (papa2 + microscape) and R (dada2)
# environments so the pipeline runs without conda env creation.
# ============================================================
FROM condaforge/miniforge3:latest

LABEL org.opencontainers.image.source="https://github.com/rec3141/microscape-nf"
LABEL org.opencontainers.image.description="Microscape amplicon sequencing pipeline"

# System dependencies
# - default-jre-headless: Nextflow needs a JVM at both build time (nextflow
#   -version below) and runtime; the miniforge3 base ships no Java.
# - build-essential + zlib1g-dev: papa2 compiles a C lib (`make libpapa2.so`)
#   during its wheel build; it needs make + gcc/g++ and the zlib dev headers.
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps curl default-jre-headless build-essential zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow and bake its framework jar into the image.
# The launcher self-downloads the framework jar on first run; on offline HPC
# compute nodes that fails. Caching under $HOME/.nextflow is not portable
# because apptainer overrides HOME at runtime (to the host home), orphaning the
# cache — so pin NXF_HOME to a fixed in-image path and bake the jar there. The
# second, network-free invocation verifies the bake so the build fails loudly
# if the offline path ever regresses.
ENV NXF_HOME=/opt/nextflow
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow -version \
    && chmod -R a+rX "$NXF_HOME" \
    && NXF_OFFLINE=true nextflow -version

# Copy environment specs
COPY envs/ /tmp/envs/

# Create Python environment (papa2 + microscape from bioconda)
# No `mamba clean` here: /opt/conda/pkgs is a BuildKit cache mount (already
# excluded from the image layer), and cleaning it fails with "Device or
# resource busy" because it's a live mount point.
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -y -p /opt/conda/envs/microscape-python \
        -f /tmp/envs/python.yml

# Create R environment (dada2, DECIPHER from bioconda)
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -y -p /opt/conda/envs/microscape-r \
        -f /tmp/envs/r.yml

# Put both envs on PATH so tools are available without conda activate
ENV PATH="/opt/conda/envs/microscape-python/bin:/opt/conda/envs/microscape-r/bin:${PATH}"

# Copy pipeline code
COPY . /pipeline/
RUN chmod +x /pipeline/bin/*.R /pipeline/entrypoint.sh

# Create docker.config that disables per-process conda (tools are on PATH)
RUN printf 'conda.enabled = false\n' > /pipeline/docker.config

# Strip caches to reduce image size
RUN find /opt/conda -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true \
    && find /opt/conda -name "*.pyc" -delete 2>/dev/null || true \
    && rm -rf /tmp/envs

# Smoke test
RUN python3 -c "import papa2; print('papa2', papa2.__version__)" \
    && python3 -c "import microscape; print('microscape', microscape.__version__)" \
    && cutadapt --version

# Writable home for non-root users (HPC)
RUN mkdir -p /home/microscape && chmod 777 /home/microscape
ENV HOME=/home/microscape

WORKDIR /data
ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
