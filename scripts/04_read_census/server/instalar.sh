#!/bin/bash
# Instala DIAMOND 2.1.10 y sra-tools en ./bin sin permisos de administrador y construye la base de la familia DCX/DCLK.
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p bin
if [ ! -x bin/diamond ]; then
  curl -sL https://github.com/bbuchfink/diamond/releases/download/v2.1.10/diamond-linux64.tar.gz | tar xz -C bin
fi
if ! ls -d bin/sratoolkit*/bin >/dev/null 2>&1; then
  curl -sL https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz | tar xz -C bin
fi
SRA=$(ls -d bin/sratoolkit*/bin | head -1)
"$SRA/vdb-config" --report-cloud-identity no >/dev/null 2>&1 || true
bin/diamond makedb --in famdb.faa -d fam --quiet
bin/diamond version; "$SRA/fastq-dump" --version | head -2; python3 --version
echo "Instalacion lista. Prueba de red con NCBI:"
"$SRA/fastq-dump" -X 10 --fasta 0 --stdout ERR3284980 | head -2
