# syntax=docker/dockerfile:1.6
FROM ghcr.io/percolator/percolator:master

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 python3-pip ca-certificates \
        libglib2.0-0 libgomp1 && \
    rm -rf /var/lib/apt/lists/* && \
    python3 -m pip install --break-system-packages --no-cache-dir \
        numpy pandas matplotlib pyIsoPEP pyprophet

COPY context /opt/context-pkg/context
ENV PYTHONPATH=/opt/context-pkg

WORKDIR /work

ENTRYPOINT ["python3", "-m", "context"]
CMD ["--help"]
