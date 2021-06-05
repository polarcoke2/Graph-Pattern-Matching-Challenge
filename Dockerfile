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

# cd does not work in Dockerfile. Use WORKDIR instead
WORKDIR /app
COPY . .
RUN mkdir build
# the result of cmake is located in the work directory thus WORKDIR /app/build
WORKDIR /app/build
RUN cmake /app
RUN make 
# change the current directory vack to /app
WORKDIR /app
CMD ["/bin/sh"]
