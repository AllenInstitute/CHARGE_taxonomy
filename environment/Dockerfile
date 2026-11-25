# hash:sha256:a954aab69542710c14c553ba11f6eb275960e62d911376f1badcc559ec970aff
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
        git=1:2.34.1-1ubuntu1.18 \
        libcurl4-openssl-dev=7.81.0-1ubuntu1.30 \
        libglpk-dev=5.0-1 \
        libgmp-dev=2:6.2.1+dfsg-3ubuntu1 \
        libhdf5-dev=1.10.7+repack-4ubuntu2 \
        liblapack-dev=3.10.0-2ubuntu1 \
        libopenblas-dev=0.3.20+ds-1 \
        libssl-dev=3.0.2-0ubuntu1.30 \
        libxml2-dev=2.9.13+dfsg-1ubuntu0.14 \
        zlib1g-dev=1:1.2.11.dfsg-2ubuntu9.2 \
    && rm -rf /var/lib/apt/lists/*

ADD "https://github.com/coder/code-server/releases/download/v4.100.3/code-server-4.100.3-linux-amd64.tar.gz" /.code-server/code-server.tar.gz

RUN cd /.code-server \
    && tar -xvf code-server.tar.gz \
    && rm code-server.tar.gz \
    && ln -s /.code-server/code-server-4.100.3-linux-amd64/bin/code-server /usr/bin/code-server \
    && mkdir -p /.vscode/extensions

RUN code-server --extensions-dir="/.vscode/extensions" --install-extension ms-python.python@2025.14.0 \
    && code-server --extensions-dir="/.vscode/extensions" --install-extension ms-toolsai.jupyter@2025.4.1 \
    && code-server --extensions-dir="/.vscode/extensions" --install-extension reditorsupport.r@2.8.6 \
    && code-server --extensions-dir="/.vscode/extensions" --install-extension saoudrizwan.claude-dev@3.35.1

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
