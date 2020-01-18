FROM nfcore/base:1.7
LABEL authors="Chris Fields" \
      description="Docker image containing all requirements for nf-core/hpcmeta pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-hpcmeta-1.0dev/bin:$PATH
