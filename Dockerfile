FROM ubuntu:18.04

LABEL maintainer="Sakkayaphab Piwluang <sakkayaphab@gmail.com>"

WORKDIR /project
COPY . .

RUN apt-get update
RUN apt-get install cmake libhts-dev libtbb-dev -y
RUN mkdir build && cd build && cmake .. -DINCLUDE_LIBRARY_PREFIX=/usr/include -DLIBRARY_LINK_PREFIX=/usr/lib/x86_64-linux-gnu
RUN cd build && make
RUN cd build && make install
RUN bolt