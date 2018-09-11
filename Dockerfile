FROM nfcore/base
LABEL
    authors="Mamana.Mbiyavanga@uct.ac.za, ayton.meintjes@uct.ac.za" \
    description="Docker image containing all requirements for h3achipimputation pipeline" \
    maintainer="Mamana Mbiyavanga <mamana.mbiyavanga@uct.ac.za>, Ayton Meintjes <ayton.meintjes@uct.ac.za>"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/h3achipimputation-1.0/bin:$PATH
