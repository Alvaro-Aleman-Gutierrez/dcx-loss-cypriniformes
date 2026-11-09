#!/bin/bash
# Repetir RELAX con variación sinónima de tasas (--srv Yes) en dclk2a y zfpl1 (K sospechosos: 23 y 8,7).
P=${PROJECT}; O=${DATA}/srv; H=${HOME}/anaconda3/envs/hyphy/bin/hyphy
cd $O; : > srv_summary.tsv
for g in dclk2a zfpl1; do for s in todas_las_perdidas solo_Danio; do
  t=$g.$s; [ -s $t.srv.json ] || $H relax --alignment $P/relax/$t.fas --tree $P/relax/$t.nwk --test Test --reference Reference \
     --models Minimal --srv Yes --output $O/$t.srv.json CPU=2 > $t.srv.log 2>&1
  python3 -I -c "import json,sys;d=json.load(open('$t.srv.json'))['test results'];print('$g','$s','srv',d['relaxation or intensification parameter'],d['p-value'],sep='\t')" >> srv_summary.tsv
done; done
echo FIN >> srv_summary.tsv
