# Convenience wrappers. `make` targets mirror exactly what CI runs, so a green
# local run predicts a green pipeline.

# Point this at a prebuilt base image to skip recompiling bcftools locally.
BASE_IMAGE ?= lrma-aou2-panel-creation-rust:latest
APP_IMAGE  ?= lrma-aou2-panel-creation-tools:latest

.PHONY: all build test fmt clippy lock docker docker-test clean

all: build test

# Debug build is faster and sufficient for the integration tests.
build:
	cargo build --workspace --locked

# Exact-match integration suite against the freshly built debug binaries.
test: build
	resources/tests/run_all.sh --bin-dir target/debug

fmt:
	cargo fmt --all -- --check

clippy:
	cargo clippy --workspace --all-targets --locked -- -D warnings

# One-time (and after any dependency change): create/refresh the committed lock.
lock:
	cargo generate-lockfile

# Build the final image containing the three release binaries. Pass a prebuilt
# base to avoid recompiling bcftools:  make docker BASE_IMAGE=...:tag
docker:
	docker build -f docker/app/Dockerfile \
		--build-arg BASE_IMAGE=$(BASE_IMAGE) \
		-t $(APP_IMAGE) .

# Run the suite *inside* the final image, proving the shipped binaries pass.
# Tests are bind-mounted, never baked in.
docker-test: docker
	docker run --rm -v "$(CURDIR)/resources/tests:/tests:ro" \
		$(APP_IMAGE) bash /tests/run_all.sh --bin-dir /usr/local/bin

clean:
	cargo clean
