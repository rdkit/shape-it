FROM continuumio/miniconda3
MAINTAINER Greg Landrum <greg.landrum@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get update && apt-get upgrade -y && apt install -y \
  curl \
  wget \
  git \
  make 

ENV LANG C

RUN conda update -n base -c defaults conda
RUN conda install gxx_linux-64 
RUN conda install -c conda-forge rdkit cmake

RUN mkdir /src
WORKDIR /src
RUN git clone https://github.com/rdkit/shape-it.git
RUN mkdir shape-it/build
WORKDIR shape-it/build
SHELL ["/bin/bash", "--login", "-c"]
RUN conda init bash
RUN conda activate base \
  && cmake -DRDKIT_INCLUDE_DIR=/opt/conda/include/rdkit -D BUILD_RDKIT_SUPPORT=ON -DBUILD_PYTHON_SUPPORT=ON ..\
  && make -j4 \
  && ctest --output-on-failure \
  && cd ../pyshapeit \
  && python -m unittest basic_test.py


