
FROM ubuntu

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y git python3-full python3-venv build-essential pipx zlib1g-dev && \
    apt-get clean

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install xgboost scikit-learn pandas

RUN git clone https://github.com/yutingsdu/TransGram && \
    cd TransGram && \
    make all release

WORKDIR /TransGram/sample_test


