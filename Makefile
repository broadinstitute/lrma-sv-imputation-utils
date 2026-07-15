# Convenience wrappers. `make` targets mirror exactly what CI runs, so a green
# local run predicts a green pipeline.
#
# The three tools are INDEPENDENT crates, each with its own committed Cargo.lock
# (reproducing the original builds). They are built separately with --locked and
# their release binaries collected into ./dist so the test harness and the
# Docker image see them in one place.

BASE_IMAGE ?= lrma-aou2-panel-creation-rust:latest
APP_IMAGE  ?= lrma-aou2-panel-creation-tools:latest

CRATES = extract-bubble-PLs pop-glimpse2 paste-vcfs

# crate -> produced binary name (default Cargo bin naming)
BIN_extract-bubble-PLs = extract-bubble-PLs
BIN_pop-glimpse2       = pop-glimpse2-joint-opt
BIN_paste-vcfs         = paste-vcfs

.PHONY: all build dist test fmt clippy lock docker docker-test clean

all: test

# Build every crate against its own committed lockfile.
build:
	@for c in $(CRATES); do \
		echo ">> building $$c"; \
		cargo build --release --locked --manifest-path resources/$$c/Cargo.toml || exit $$?; \
	done

# Collect the release binaries into ./dist (one directory, like the image).
dist: build
	@mkdir -p dist
	@cp resources/extract-bubble-PLs/target/release/extract-bubble-PLs dist/
	@cp resources/pop-glimpse2/target/release/pop-glimpse2-joint-opt   dist/
	@cp resources/paste-vcfs/target/release/paste-vcfs                 dist/
	@echo ">> dist/: $$(ls dist)"

# Exact-match integration suite against the collected binaries.
test: dist
	bash resources/tests/run_all.sh --bin-dir dist

fmt:
	@for c in $(CRATES); do cargo fmt --manifest-path resources/$$c/Cargo.toml -- --check || exit $$?; done

clippy:
	@for c in $(CRATES); do \
		cargo clippy --locked --all-targets --manifest-path resources/$$c/Cargo.toml -- -D warnings || exit $$?; \
	done

# Regenerate the per-crate lockfiles (rarely needed: the committed locks are the
# reproducibility contract). Only run this when intentionally updating deps.
lock:
	@for c in $(CRATES); do cargo generate-lockfile --manifest-path resources/$$c/Cargo.toml; done

docker:
	docker build -f docker/app/Dockerfile \
		--build-arg BASE_IMAGE=$(BASE_IMAGE) \
		-t $(APP_IMAGE) .

# Run the suite inside the final image; tests bind-mounted, never baked in.
docker-test: docker
	docker run --rm -v "$(CURDIR)/resources/tests:/tests:ro" \
		$(APP_IMAGE) bash /tests/run_all.sh --bin-dir /usr/local/bin

clean:
	@for c in $(CRATES); do cargo clean --manifest-path resources/$$c/Cargo.toml; done
	rm -rf dist
