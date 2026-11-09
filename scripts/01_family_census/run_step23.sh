#!/bin/bash
# Orquestador pasos 2-3: descarga genomas de referencia y ejecuta step23_pipeline.py
P=${PROJECT}
echo "[$(date)] descarga list3" >> $P/step23.log
bash $P/genomes/dl.sh $P/genomes/list3.txt >> $P/step23.log 2>&1
echo "[$(date)] pipeline" >> $P/step23.log
python3 -I $P/scripts/step23_pipeline.py >> $P/step23.log 2>&1
echo "[$(date)] FIN (exit $?)" >> $P/step23.log
