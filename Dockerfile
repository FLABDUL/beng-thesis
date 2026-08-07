FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
       build-essential \
       ca-certificates \
       cmake \
       libeigen3-dev \
       libpcl-dev \
       ninja-build \
       zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
COPY . .

RUN cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build --parallel \
    && cmake --install build --prefix /opt/beng-thesis

ENV PATH="/opt/beng-thesis/bin:${PATH}"

CMD ["compute_ma", "--help"]
