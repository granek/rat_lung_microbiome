## Derived from https://raw.githubusercontent.com/jupyter/docker-stacks/master/base-notebook/Dockerfile

# Distributed under the terms of the Modified BSD License.

FROM rocker/rstudio

MAINTAINER Josh Granek

USER root


# Configure environment
ENV CONDA_DIR /opt/conda
ENV PATH $CONDA_DIR/bin:$PATH
ENV SHELL /bin/bash
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV RSTUDIO_USER rstudio

# Create jovyan user with UID=1000 and in the 'users' group
RUN mkdir -p $CONDA_DIR && \
    chown $RSTUDIO_USER $CONDA_DIR

USER $RSTUDIO_USER



# Install conda as ?????
RUN cd /tmp && \
    mkdir -p $CONDA_DIR && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.1.11-Linux-x86_64.sh && \
    echo "efd6a9362fc6b4085f599a881d20e57de628da8c1a898c08ec82874f3bad41bf *Miniconda3-4.1.11-Linux-x86_64.sh" | sha256sum -c - && \
    /bin/bash Miniconda3-4.1.11-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm Miniconda3-4.1.11-Linux-x86_64.sh && \
    $CONDA_DIR/bin/conda install --quiet --yes conda==4.1.11 && \
    $CONDA_DIR/bin/conda config --system --add channels conda-forge && \
    $CONDA_DIR/bin/conda config --system --set auto_update_conda false && \
    conda clean -tipsy

# Install qiime1 notebook as 
RUN conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda && \
    conda clean -tipsy

USER root
CMD ["/init"]