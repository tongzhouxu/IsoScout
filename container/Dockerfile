FROM mambaorg/micromamba:2.0.5

COPY --chown=$MAMBA_USER:$MAMBA_USER container/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/isoscout
ENV PATH="/opt/isoscout/scripts:${PATH}"
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python3", "/opt/isoscout/scripts/screen_isolate.py"]
