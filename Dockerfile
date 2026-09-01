FROM nbeck/mfannot

RUN apt-get update && apt-get install -y \
    aragorn \
    mafft \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

COPY src/ /opt/pipeline/src/
WORKDIR /opt/pipeline/src

ENTRYPOINT ["python3", "run.py"]