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
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps curl \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow -version

# Copy environment specs
COPY envs/ /tmp/envs/

# Create Python environment (papa2 + microscape from bioconda)
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -y -p /opt/conda/envs/microscape-python \
        -f /tmp/envs/python.yml \
    && mamba clean -afy

# Create R environment (dada2, DECIPHER from bioconda)
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -y -p /opt/conda/envs/microscape-r \
        -f /tmp/envs/r.yml \
    && mamba clean -afy

# Put both envs on PATH so tools are available without conda activate
ENV PATH="/opt/conda/envs/microscape-python/bin:/opt/conda/envs/microscape-r/bin:${PATH}"

# Provenance, baked at build time.
#
# There is no git repository inside the image, so `workflow.commitId` is null for every
# run launched from the .sif — which is how this pipeline actually runs. Three different
# main.nf versions were once in play at once (the one that produced a dataset, the one
# shipped in the container, and the one in git) with no way to tell them apart. CI
# already tags the image with the SHA; this puts it somewhere the pipeline can read.
ARG GIT_SHA=unknown
ARG GIT_REF=unknown
ARG BUILD_DATE=unknown
ENV MICROSCAPE_GIT_SHA=${GIT_SHA} \
    MICROSCAPE_GIT_REF=${GIT_REF} \
    MICROSCAPE_BUILD_DATE=${BUILD_DATE}
LABEL org.opencontainers.image.revision="${GIT_SHA}" \
      org.opencontainers.image.version="${GIT_REF}" \
      org.opencontainers.image.created="${BUILD_DATE}"

# Copy pipeline code
COPY . /pipeline/
RUN chmod +x /pipeline/bin/*.R /pipeline/bin/*.sh /pipeline/entrypoint.sh

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
