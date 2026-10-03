# syntax=docker/dockerfile:1
FROM nvidia/cuda:12.6.3-devel-ubuntu24.04 AS toolchain
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake ninja-build g++ libspdlog-dev zlib1g-dev catch2 python3 \
    && rm -rf /var/lib/apt/lists/*
FROM toolchain AS compile
WORKDIR /src
COPY . .
ARG CUDA_ARCHITECTURES="70;75;80;86;89;90"
ARG BUILD_JOBS=2
RUN cmake -S . -B build-container -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CUDA_ARCHITECTURES="${CUDA_ARCHITECTURES}" \
    && cmake --build build-container -j "${BUILD_JOBS}"

FROM compile AS build
RUN ctest --test-dir build-container --output-on-failure \
    && build-container/src/stomata --quickstart --cpu-only \
    && cmake --install build-container --prefix /opt/stomata

FROM nvidia/cuda:12.6.3-runtime-ubuntu24.04 AS runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
    libspdlog1.12 zlib1g libstdc++6 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/stomata/bin/stomata /usr/local/bin/stomata
WORKDIR /work
ENTRYPOINT ["stomata"]
CMD ["--help"]
