FROM centos:centos7

RUN yum update -y
RUN yum install -y wget make
RUN mkdir -p /tmp/cmake && \
  wget 'https://cmake.org/files/v3.20/cmake-3.20.2-linux-x86_64.sh' && \
  bash cmake-3.20.2-linux-x86_64.sh --prefix=/usr/local --exclude-subdir && \
  rm -rf /tmp/cmake
RUN yum install -y centos-release-scl
RUN yum install -y devtoolset-7-gcc*
ENV PATH="/opt/rh/devtoolset-7/root/usr/bin:$PATH"
RUN source scl_source enable devtoolset-7
RUN yum clean all

WORKDIR /app
COPY . .
RUN mkdir build
WORKDIR /app/build
RUN cmake /app
RUN make 
CMD ["/bin/sh"]
