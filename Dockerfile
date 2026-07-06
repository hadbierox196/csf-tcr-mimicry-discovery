# syntax=docker/dockerfile:1
#
# Two build targets from one file:
#   docker build --target scoring .              (default: CPU-only, lightweight)
#   docker build --target structure-prediction .  (CUDA-enabled; heavier)
#
# Model weights (multi-GB) are intentionally NOT baked into either image --
# they're DVC-tracked (see models/weights.dvc) and pulled or mounted at
# container start (see docker-compose.yml). Baking multi-GB weights into an
# image layer bloats every rebuild and duplicates what DVC already versions.
#
# For HPC nodes that don't allow rootful Docker, convert the published image
# with containers/singularity/build_sif.sh instead of using this file directly.

########################################################################
# Target: scoring -- CPU-only, default. io + scoring + CLI.
########################################################################
FROM python:3.11-slim AS scoring

LABEL org.opencontainers.image.title="mimicry-discovery (scoring)"
LABEL org.opencontainers.image.description="CPU-only image: io, scoring, and the CLI. No GPU or structure-prediction dependencies."

RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/mimicry-discovery
COPY pyproject.toml README.md ./
COPY mimicry_discovery/ ./mimicry_discovery/
RUN pip install --no-cache-dir .

RUN useradd --create-home --uid 1000 mimicry
USER mimicry

ENTRYPOINT ["mimicry-discovery"]
CMD ["--help"]

########################################################################
# Target: structure-prediction -- CUDA-enabled.
#
# TCRdock is now wired in (roadmap Sprint 2, see
# mimicry_discovery/structure/tcrdock_adapter.py) -- this stage clones
# the real TCRdock repo. What's still missing, deliberately, because it
# needs setup this Dockerfile doesn't own or fully know the shape of:
#   - AlphaFold itself (TCRdock runs AlphaFold under the hood) -- its
#     genetic-database + parameter setup is extensive and versioned
#     independently; follow TCRdock's own README for the current steps
#     rather than duplicating them here where they'd drift out of sync.
#   - ESMFold / AlphaFold-Multimer's own dependencies (torch, etc.) --
#     still scaffolded (roadmap Sprint 3), so not installed yet for the
#     same reason as before: no benefit to a multi-GB wheel for code
#     that can't call it.
########################################################################
FROM nvidia/cuda:12.6.1-cudnn-runtime-ubuntu24.04 AS structure-prediction

LABEL org.opencontainers.image.title="mimicry-discovery (structure-prediction)"
LABEL org.opencontainers.image.description="CUDA image with TCRdock cloned; AlphaFold setup and ESMFold/AF-Multimer deps still required/pending -- see the comment above."

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip python3-venv git wget \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/bin/python3 /usr/bin/python

WORKDIR /opt/mimicry-discovery
COPY pyproject.toml README.md ./
COPY mimicry_discovery/ ./mimicry_discovery/
RUN pip install --no-cache-dir --break-system-packages .

# TCRdockPredictor defaults to /opt/TCRdock (see tcrdock_adapter.py) --
# clone at a pinned commit, not a moving branch, so a rebuild months from
# now doesn't silently pick up upstream changes.
# TODO(containers): pin an actual commit SHA once Sprint 2 is confirmed
# working end to end against a real AlphaFold install; "main" below is a
# placeholder for the first real build, not a long-term pin.
RUN git clone --depth 1 https://github.com/phbradley/TCRdock.git /opt/TCRdock

# TODO(containers): follow TCRdock's own README for AlphaFold setup
# (genetic databases + model parameters) -- out of scope to replicate
# here since it's large, versioned independently, and best kept as the
# upstream source of truth rather than a copy that can drift.

# TODO(containers): once ESMFold/AlphaFold-Multimer adapters are
# implemented (roadmap Sprint 3), add their deps here, e.g.:
#   RUN pip install --no-cache-dir --break-system-packages \
#         torch --index-url https://download.pytorch.org/whl/cu124
#   RUN pip install --no-cache-dir --break-system-packages fair-esm

RUN useradd --create-home --uid 1000 mimicry
USER mimicry

ENTRYPOINT ["mimicry-discovery"]
CMD ["--help"]
