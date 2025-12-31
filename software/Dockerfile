FROM python:3.10.19

ENV GIT_PYTHON_REFRESH=quiet

RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    build-essential \
    libgomp1 \
    gnupg \
    ca-certificates \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get install ca-certificates gnupg
RUN gpg --homedir /tmp --no-default-keyring --keyring gnupg-ring:/usr/share/keyrings/mono-official-archive-keyring.gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN chmod +r /usr/share/keyrings/mono-official-archive-keyring.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/mono-official-archive-keyring.gpg] https://download.mono-project.com/repo/ubuntu stable-focal main" | tee /etc/apt/sources.list.d/mono-official-stable.list
RUN apt-get update && apt-get install -y mono-devel

RUN wget https://github.com/MannLabs/alphapept/blob/master/alphapept/ext/bruker/FF/linux64/libtbb.so.2 \
    && mv libtbb.so.2 /usr/lib/

WORKDIR /app

RUN pip install uv -i https://mirrors.aliyun.com/pypi/simple/

RUN uv venv /app/.venv

COPY . /app/

RUN uv pip install -e . --refresh -i https://mirrors.aliyun.com/pypi/simple/

ENV PATH="/app/.venv/bin:$PATH"

ENTRYPOINT ["python", "-m", "dda_bert.cli"]
CMD ["--help"]