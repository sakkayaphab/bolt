FROM ubuntu:18.04

LABEL maintainer="Sakkayaphab Piwluang <sakkayaphab@gmail.com>"

WORKDIR /project
COPY . .

RUN apt-get update
