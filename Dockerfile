FROM nfcore/base
LABEL authors="aytonm@gmail.com;mamana.mbiyavanga@uct.ac.za" \
      description="for now to keep nf-core happy"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-chipimputation-1.0dev/bin:$PATH
