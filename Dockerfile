FROM ubuntu:18.04

LABEL maintainer="Sakkayaphab Piwluang <sakkayaphab@gmail.com>"

WORKDIR /project
COPY . .

RUN apt-get -y update
RUN apt-get -y install build-essential
RUN apt-get -y install manpages-dev
RUN gcc --version
RUN apt-get -y install cmake libhts-dev libtbb-dev -y
RUN mkdir build && cd build && cmake .. -DINCLUDE_LIBRARY_PREFIX=/usr/include -DLIBRARY_LINK_PREFIX=/usr/lib/
RUN cd build && make
RUN cd build && make install

CMD ["bolt"]
