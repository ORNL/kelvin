FROM fedora:30

RUN yum update -y
RUN yum install -y git cmake make gcc-c++ boost-devel libstdc++-devel doxygen

RUN git clone https://github.com/jayjaybillings/parsers
RUN mkdir parsers-build
RUN cd parsers-build \
  && cmake ../parsers -DCMAKE_BUILD_TYPE=Debug -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_ECLIPSE_VERSION=4.5 \
  && make \
  && make test \
  && make doc \
  && make install
  
RUN git clone https://github.com/mfem/mfem.git
RUN mkdir mfem/build 
RUN cd mfem/build \
  && cmake .. \
  && make install
  
RUN git clone https://github.com/ORNL/kelvin.git
RUN mkdir kelvin/build
RUN cd kelvin/build \
  && cmake ../ -DPARSERS_DIR=/parsers -DMFEM_DIR=/usr/local \
  && make \
  && make test
