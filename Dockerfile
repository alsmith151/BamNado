# Multi-stage build for minimal image size
FROM rust:latest AS builder

WORKDIR /build

# Copy workspace files
COPY Cargo.toml Cargo.lock ./
COPY bamnado ./bamnado
COPY bamnado-python ./bamnado-python

# Build the binary in release mode
RUN cargo build --release --package bamnado

# Final stage: use ubuntu for good glibc compatibility
FROM ubuntu:24.04

# Remove unnecessary packages to reduce image size
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Copy the binary from builder
COPY --from=builder /build/target/release/bamnado /usr/local/bin/bamnado

ENTRYPOINT ["/usr/local/bin/bamnado"]
CMD ["--help"]
