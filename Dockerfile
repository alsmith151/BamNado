# Multi-stage build for minimal image size and faster builds

# Stage 1: Planner - generate recipe for cargo-chef
FROM rust:1.84 AS planner
WORKDIR /build
RUN cargo install cargo-chef
COPY Cargo.toml ./
COPY bamnado ./bamnado
COPY bamnado-python ./bamnado-python
RUN cargo chef prepare --recipe-path recipe.json

# Stage 2: Cacher - build and cache dependencies
FROM rust:1.84 AS cacher
WORKDIR /build
RUN cargo install cargo-chef
COPY --from=planner /build/recipe.json recipe.json
RUN cargo chef cook --release --recipe-path recipe.json

# Stage 3: Builder - build the application
FROM rust:1.84 AS builder
WORKDIR /build
COPY Cargo.toml ./
COPY bamnado ./bamnado
COPY bamnado-python ./bamnado-python
COPY --from=cacher /build/target target
RUN cargo build --release --package bamnado

# Stage 4: Runtime - minimal final image
FROM debian:bookworm-slim

# Install minimal dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Copy the binary from builder
COPY --from=builder /build/target/release/bamnado /usr/local/bin/bamnado

ENTRYPOINT ["/usr/local/bin/bamnado"]
CMD ["--help"]
