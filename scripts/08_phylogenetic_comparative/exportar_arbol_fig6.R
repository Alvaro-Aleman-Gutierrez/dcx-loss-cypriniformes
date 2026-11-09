#!/usr/bin/env Rscript
# Exporta el arbol de Sudasinghe et al. 2026 podado a las puntas del censo (mismas reglas que analisis_filogenetico.R)
# para la Fig. 6 (figuras.py fig6): integracion/arbol_censo.tre + integracion/arbol_censo_puntas.tsv
suppressMessages({ .libPaths(c("~/R/libs", .libPaths())); library(ape) })
D <- "${DATA}"; O <- file.path(D, "integracion")
tr <- read.nexus(file.path(D, "sudasinghe2026/Sudasinghe-etal.-BEAST-time.calibrated.phylogeny-335taxa.tre"))
if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
T <- read.delim(file.path(O, "tabla_dcx_arbol.tsv"), stringsAsFactors = FALSE)
T <- T[!is.na(T$punta) & T$punta != "" & T$estado %in% c("presente", "ausente"), ]
T <- T[order(T$punta, T$fuente != "lecturas"), ]; T <- T[!duplicated(T$punta), ]
T <- T[!grepl("^Characiformes|^Siluriformes|^Gymnotiformes", T$punta), ]
tr <- ladderize(keep.tip(tr, intersect(tr$tip.label, T$punta)))
write.tree(tr, file.path(O, "arbol_censo.tre"))
write.table(T[match(tr$tip.label, T$punta), c("punta", "especie", "fuente", "estado", "familia", "miniatura")],
            file.path(O, "arbol_censo_puntas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(Ntip(tr), "puntas\n")
