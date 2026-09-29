# hash:sha256:65c7a47e2754c6bb848be7311b62c806b2eca4c539aaebf085f741064b52e75e
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/py-r:python3.10.12-R4.2.3-IRkernel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

ENV CRAN=https://packagemanager.posit.co/cran/__linux__/jammy/latest
ENV LANG=en_US.UTF-8
ENV TZ=Etc/UTC

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential=12.9ubuntu3 \
        gfortran=4:11.2.0-1ubuntu1 \
        git=1:2.34.1-1ubuntu1.17 \
        libcurl4-openssl-dev=7.81.0-1ubuntu1.29 \
        libglpk-dev=5.0-1 \
        libgmp-dev=2:6.2.1+dfsg-3ubuntu1 \
        libhdf5-dev=1.10.7+repack-4ubuntu2 \
        liblapack-dev=3.10.0-2ubuntu1 \
        libopenblas-dev=0.3.20+ds-1 \
        libssl-dev=3.0.2-0ubuntu1.29 \
        libxml2-dev=2.9.13+dfsg-1ubuntu0.13 \
        zlib1g-dev=1:1.2.11.dfsg-2ubuntu9.2 \
    && rm -rf /var/lib/apt/lists/*

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
