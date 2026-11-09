#!/bin/bash
# Arbol de la Figura 1 (matriz extendida), en el orden metodologicamente correcto:
#  1) SELECCION DE MODELO: particiones de partida gen x posicion de codon (49 x 3 = 147);
#     ModelFinder + MERGE (-m MF+MERGE, relaxed clustering -rcluster 10) elige el modelo de cada
#     particion y fusiona las de dinamica similar. Solo seleccion, sin busqueda de arbol final.
#  2) 10 BUSQUEDAS ML INDEPENDIENTES (--runs 10) con el esquema elegido (best_scheme.nex),
#     1000 UFBoot + 1000 SH-aLRT, grupo externo Astyanax + Ictalurus.
#  3) Convergencia: logL de cada corrida, RF entre las 10, y RF frente al arbol GTR+F+G4 (ml.treefile).
# (2026-11-05: sustituye a arbol_ext_runs10.sh, que hacia las 10 corridas con GTR+F+G4 sin seleccionar.)
set -uo pipefail
O=${DATA}/arbol_ext
T=${THREADS:-10}
python3 - "$O" <<'EOF'
import re, sys
O = sys.argv[1]
cs = re.findall(r"charset (\S+) = (\d+)-(\d+);", open(f"{O}/partitions.nex").read())
with open(f"{O}/partitions_codpos.nex", "w") as fo:
    fo.write("#nexus\nbegin sets;\n")
    for g, a, b in cs:
        a, b = int(a), int(b)
        assert (b - a + 1) % 3 == 0, f"{g}: longitud no multiplo de 3"
        for k in range(3):
            fo.write(f"  charset {g}_p{k+1} = {a+k}-{b}\\3;\n")
    fo.write("end;\n")
print(f"particiones gen x posicion: {3*len(cs)}")
EOF
echo "[$(date)] 1) ModelFinder + MERGE con $T hilos"
nice -n 10 iqtree2 -s $O/super.fas -p $O/partitions_codpos.nex -m MF+MERGE -rcluster 10 \
  -T $T --prefix $O/mf -redo
[ -s $O/mf.best_scheme.nex ] || { echo "FALLO: no hay mf.best_scheme.nex"; exit 1; }
echo "== esquema elegido: $(grep -c 'charset' $O/mf.best_scheme.nex) particiones tras fusionar"
grep -E "charpartition" -A200 $O/mf.best_scheme.nex | tr ',' '\n' | sed 's/:.*//' | sort | uniq -c | sort -rn | head -15
echo "[$(date)] 2) --runs 10 con el esquema elegido"
nice -n 10 iqtree2 -s $O/super.fas -p $O/mf.best_scheme.nex -B 1000 -alrt 1000 \
  --runs 10 -T $T --prefix $O/mf10 -redo -o Astyanax_mexicanus,Ictalurus_furcatus
echo "== logL por corrida"
grep -E "BEST SCORE FOUND|START RUN NUMBER" $O/mf10.log
echo "== RF entre las 10 corridas (0 = misma topologia)"
iqtree2 -rf_all $O/mf10.runtrees --prefix $O/mf10_rf -redo >/dev/null 2>&1
cat $O/mf10_rf.rfdist
echo "== RF: GTR+F+G4 (ml) vs mejor de las 10 con modelo seleccionado (mf10)"
cat $O/ml.treefile $O/mf10.treefile > $O/gtr_vs_mf.trees
iqtree2 -rf_all $O/gtr_vs_mf.trees --prefix $O/gtr_vs_mf -redo >/dev/null 2>&1
cat $O/gtr_vs_mf.rfdist
echo "FIN mfp"
