FROM mambaorg/micromamba:1.5.8

USER root
RUN micromamba install -y -n base \
    -c conda-forge -c bioconda \
    python=3.11 r-base gapseq=1.4.0 \
    r-data.table r-jsonlite r-stringr r-optparse r-sybil \
    && micromamba clean --all --yes

WORKDIR /app
COPY . /app
RUN python -m pip install --no-cache-dir .
ENV PATH=/opt/conda/bin:/opt/micromamba/bin:$PATH
ENTRYPOINT ["gemagent"]