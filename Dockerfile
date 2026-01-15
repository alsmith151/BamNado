# Multi-stage build for minimal image size

# Stage 1: Builder
FROM rust:nightly AS builder
WORKDIR /build

# Copy workspace files
COPY Cargo.toml ./
COPY bamnado ./bamnado
COPY bamnado-python ./bamnado-python

# Build the binary in release mode
RUN cargo build --release --package bamnado

# Stage 2: Runtime - minimal final image
FROM debian:bookworm-slim

# Install minimal dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Copy the binary from builder
COPY --from=builder /build/target/release/bamnado /usr/local/bin/bamnado

ENTRYPOINT ["/usr/local/bin/bamnado"]
CMD ["--help"]
