#!/usr/bin/env Rscript
# Q1 del articulo: cuantas veces y cuando se perdio DCX, y si se asocia con la miniaturizacion.
# Arbol: Sudasinghe et al. 2026 (PNAS), BEAST, 335 taxones, Ma (Zenodo 10.5281/zenodo.23008328, CC-BY).
# Datos: integracion/tabla_dcx_arbol.tsv (scripts/integrar_censo.py). Solo Cypriniformes; estados 'ambiguo' y
# 'cobertura_baja' se excluyen.
# 1) Dollo: DCX presente en el ancestro; cada clado maximo con todas las puntas ausentes = una perdida.
#    Edad: intervalo de la rama troncal del clado (edad del nodo padre -> edad del nodo del clado; puntas = 0).
# 2) Mapeo estocastico (phytools::make.simmap, modelo irreversible presente->ausente, 100 simulaciones): numero de perdidas.
# 3) Pagel (phytools::fitPagel): DCX ausente vs (a) miniatura progenetica, (b) progenetica o enana proporcionada.
suppressMessages({ .libPaths(c("~/R/libs", .libPaths())); library(ape); library(phytools) })
set.seed(20261107)
D <- "${DATA}"; O <- file.path(D, "integracion")
tr <- read.nexus(file.path(D, "sudasinghe2026/Sudasinghe-etal.-BEAST-time.calibrated.phylogeny-335taxa.tre"))
if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
T <- read.delim(file.path(O, "tabla_dcx_arbol.tsv"), stringsAsFactors = FALSE)
T <- T[!is.na(T$punta) & T$punta != "" & T$estado %in% c("presente", "ausente"), ]
# una fila por punta (prioridad: lecturas sobre genoma; no hay conflictos, ver integrar_censo.py)
T <- T[order(T$punta, T$fuente != "lecturas"), ]; T <- T[!duplicated(T$punta), ]
T <- T[!grepl("^Characiformes|^Siluriformes|^Gymnotiformes", T$punta), ]
tr <- keep.tip(tr, intersect(tr$tip.label, T$punta))
x <- setNames(ifelse(T$estado == "ausente", "ausente", "presente"), T$punta)[tr$tip.label]
cat(sprintf("Puntas analizadas: %d (presente %d, ausente %d)\n", Ntip(tr), sum(x == "presente"), sum(x == "ausente")))

## 1) Dollo + edades
H <- max(node.depth.edgelength(tr)); edad <- H - node.depth.edgelength(tr)     # edad de cada nodo (puntas ~ 0)
perd <- c()
for (n in c(1:Ntip(tr), (Ntip(tr) + 2):(Ntip(tr) + tr$Nnode))) {
  hijos <- if (n <= Ntip(tr)) tr$tip.label[n] else extract.clade(tr, n)$tip.label
  if (all(x[hijos] == "ausente")) {
    padre <- tr$edge[tr$edge[, 2] == n, 1]
    hp <- extract.clade(tr, padre)$tip.label
    if (!all(x[hp] == "ausente")) perd <- rbind(perd, data.frame(nodo = n, n_puntas = length(hijos),
      edad_max_Ma = round(edad[padre], 1), edad_min_Ma = round(edad[n], 1),
      puntas = paste(sub("^[^_]+_[^_]+_", "", hijos), collapse = ", ")))
  }
}
cat(sprintf("\n== Dollo: %d perdidas independientes de DCX\n", nrow(perd))); write.table(perd, file.path(O, "perdidas_dollo.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pv <- perd; pv$puntas <- substr(pv$puntas, 1, 90); print(pv[, -1], row.names = FALSE)

## 2) Mapeo estocastico, modelo irreversible
Q <- matrix(c(0, 1, 0, 0), 2, 2, dimnames = list(c("ausente", "presente"), c("ausente", "presente")))  # fila=desde, col=hacia: solo presente->ausente
sm <- make.simmap(tr, x, model = Q, nsim = 100, pi = c(0, 1), message = FALSE)
nloss <- sapply(sm, function(s) sum(sapply(s$maps, function(m) length(m) > 1 && names(m)[1] == "presente")))
cat(sprintf("\n== Mapeo estocastico (irreversible, 100 sim.): perdidas mediana %g (rango %g-%g)\n",
            median(nloss), min(nloss), max(nloss)))

## 3) Pagel
mini <- setNames(T$miniatura, T$punta)[tr$tip.label]
for (cod in c("progenetica", "progenetica_o_enana")) {
  y <- if (cod == "progenetica") ifelse(mini == "progenetic-miniature", "mini", "no") else
         ifelse(mini %in% c("progenetic-miniature", "proportioned-dwarf"), "mini", "no")
  ok <- !is.na(y); tt <- keep.tip(tr, tr$tip.label[ok])
  fp <- tryCatch(fitPagel(tt, x[tt$tip.label], setNames(y[ok], tr$tip.label[ok])[tt$tip.label]), error = function(e) e)
  if (inherits(fp, "error")) { cat("\nPagel", cod, "fallo:", conditionMessage(fp), "\n"); next }
  tab <- table(DCX = x[tt$tip.label], mini = y[ok])
  cat(sprintf("\n== Pagel (%s): LR=%.2f, p=%.3g; logL indep=%.2f dep=%.2f\n", cod, fp$lik.ratio, fp$P,
              fp$independent.logL, fp$dependent.logL)); print(tab)
}
cat("\nFIN\n")
