# Installation

Stomata targets Linux x86_64. Native Windows builds are not supported; WSL2
with NVIDIA CUDA support is an option. CPU-only builds need a C++17 compiler,
CMake 3.20+, spdlog and zlib. Tests additionally need Catch2 3.

See the root README for Ubuntu 24.04 packages and a CPU build. GPU builds
require CUDA Toolkit 11.8+ and a compatible NVIDIA driver; `--cpu-only` is a
runtime switch, distinct from `-DSTOMATA_ENABLE_CUDA=OFF` at build time.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure
build/src/stomata --quickstart
cmake --install build --prefix "$HOME/.local"
```

Default CUDA targets are `70;75;80;86;89;90` (Volta through Hopper). A compiled
target is not a claim of hardware testing on every model. Override with
`-DCMAKE_CUDA_ARCHITECTURES=75` for a smaller GTX 1660/Turing build, for example.
Newer architectures/toolkits need a build and smoke test on that hardware.

## Testing

Normal tests use synthetic data and need no reference download:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSTOMATA_ENABLE_CUDA=OFF
cmake --build build -j2
ctest --test-dir build --output-on-failure
build/src/stomata --quickstart --cpu-only
```

Catch2 3 is required when `BUILD_TESTING=ON` (the default). A missing test
dependency is an error rather than a successful build with no tests. Use
`-DBUILD_TESTING=OFF` only for an installation that deliberately omits tests.

For CUDA, omit `-DSTOMATA_ENABLE_CUDA=OFF`, run CTest on a GPU host, and run
`build/src/stomata --quickstart`. Confirm the smoke-test rows say `[GPU]`;
a CPU fallback is not a GPU validation result. The independent oracle and fuzz
checks are described in [validation/README.md](../validation/README.md).

## Container

The Dockerfile builds against CUDA 12.6.3 / Ubuntu 24.04 and runs the C++ tests
and CPU smoke test during build. Building an image does not exercise a GPU.

```bash
docker build -t stomata:local .
docker run --rm stomata:local --quickstart --cpu-only
docker run --gpus all --rm stomata:local --quickstart
docker run --gpus all --rm -v "$PWD:/work" stomata:local --genome /work/hg38.fa.st --spacer-file /work/guides.tsv --threshold 3 --pam NGG --output /work/hits.tsv
```

Docker GPU use requires NVIDIA Container Toolkit on the host. The final image
contains the CLI and runtime dependencies; the larger `build` stage contains
tests and evaluation sources. To test on a GPU:

```bash
docker build --target build -t stomata-build:local .
docker run --gpus all --rm stomata-build:local ctest --test-dir /src/build-container --output-on-failure
```

For Apptainer, convert the locally built image:

```bash
docker save stomata:local -o stomata.tar
apptainer build stomata.sif docker-archive://stomata.tar
apptainer run --nv stomata.sif --quickstart
```

No public image tag is claimed until a release image has actually been pushed.
Pin a release digest when adopting a published image in a workflow.
