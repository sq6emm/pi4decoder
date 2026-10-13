# Build/test image for pi4decoder. The decoder itself needs only a C99
# compiler and libm; python3 is here for the benchmark/plot scripts only.
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc libc6-dev make time python3 valgrind \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /work
