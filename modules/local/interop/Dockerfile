# Illumina InterOp Python bindings plus the plotting stack needed by
# collect_illumina_metrics.py. CPU-only.
#
# The colleague's repo installs these with `uv sync` from a pyproject; here they
# are baked into an image so tasks stay offline-safe on the cluster.
#
# TODO: this image is not published yet. Either add a build workflow to this
# repo, or move the module to EIT-GBI/nf-mod-interop (which gets build.yml for
# free from nf-mod-template) and pin ghcr.io/eit-gbi/nf-mod-interop instead.
FROM python:3.12-slim

ARG INTEROP_VERSION=1.3.2
ARG NUMPY_VERSION=2.1.3
ARG MATPLOTLIB_VERSION=3.9.2

RUN pip install --no-cache-dir \
        "interop==${INTEROP_VERSION}" \
        "numpy==${NUMPY_VERSION}" \
        "matplotlib==${MATPLOTLIB_VERSION}"

# Pre-build the matplotlib font cache so the first task does not pay for it.
RUN MPLBACKEND=Agg python -c "import matplotlib.pyplot"
