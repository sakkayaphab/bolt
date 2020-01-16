FROM ubuntu:18.04

LABEL maintainer="Sakkayaphab Piwluang <sakkayaphab@gmail.com>"

WORKDIR /project
COPY . .

RUN apt-get -y update
RUN apt-get -y install build-essential
RUN apt-get -y install manpages-dev
RUN gcc --version
RUN apt-get -y install zlib1g-dev libncurses5-dev
RUN apt-get -y install cmake libhts-dev libtbb-dev bzip2
RUN cd / && ls -al
RUN cd /project/ && ls -al
RUN cd /usr/lib/x86_64-linux-gnu && ls -al
RUN cd /usr/local && ls -al
RUN mkdir build && cd build && cmake .. -DINSTALL_BIN_PREFIX=${PWD} -DINCLUDE_LIBRARY_PREFIX=/usr/include -DLIBRARY_LINK_PREFIX=/usr/lib/x86_64-linux-gnu
RUN cd build && make
RUN cd build && make install

CMD ["/project/build/bolt"]