# Preserve the published Mashpit/Python environment byte-for-byte.
FROM ghcr.io/tongzhouxu/isoscout@sha256:0ee316333279d2fbb28ba78605aa25ae86b0cdf9321cac15dd78f2c349ad9cf9

LABEL org.opencontainers.image.source="https://github.com/tongzhouxu/IsoScout"

COPY --chown=$MAMBA_USER:$MAMBA_USER container/alignment-linux-64.lock /tmp/alignment-linux-64.lock
RUN micromamba create -y -p /opt/conda/envs/assembly-tools -f /tmp/alignment-linux-64.lock && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/isoscout
ENV PATH="/opt/conda/envs/assembly-tools/bin:/opt/isoscout/scripts:${PATH}"
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python3", "/opt/isoscout/scripts/screen_isolate.py"]
