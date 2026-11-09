#!/usr/bin/env Rscript
# Robustez de Q1 frente al arbol: mismos analisis que analisis_filogenetico.R (Dollo, mapeo estocastico irreversible,
# Pagel) sobre arboles independientes del de Sudasinghe et al. 2026:
#  (A) Fish Tree of Life, version completa (Rabosky et al. 2018; actinopt_full.trees, 20 de las 100 replicas): especies sin
#      datos geneticos ubicadas por taxonomia. Emparejamiento por nombre exacto; si no existe (p. ej. especies no descritas
#      'sp' o D. cerebrum, descrita en 2021) se usa un congenere aun no usado (orden alfabetico) y se registra.
#  (B) Nuestro arbol ML de 24 genomas (arbol_ext/mf10.treefile), enraizado con Astyanax + Ictalurus, grupo externo eliminado
#      y hecho ultrametrico con chronos() solo para Pagel.
# Estados y miniaturizacion: integracion/tabla_dcx_arbol.tsv (una fila por especie; ambiguos excluidos).
suppressMessages({ .libPaths(c("~/R/libs", .libPaths())); library(ape); library(phytools) })
set.seed(20261107)
D <- "${DATA}"; O <- file.path(D, "integracion")
T <- read.delim(file.path(O, "tabla_dcx_arbol.tsv"), stringsAsFactors = FALSE)
T <- T[T$estado %in% c("presente", "ausente") & !T$especie %in% c("Astyanax_mexicanus", "Ictalurus_furcatus"), ]
T$especie <- sub("_sp\\..*$|_sp$", "_sp", T$especie)
T <- T[order(T$especie, T$fuente != "lecturas"), ]
# miniaturizacion por especie; si falta, la del congenere en la tabla (mismo genero)
for (i in which(is.na(T$miniatura))) {
  g <- sub("_.*", "", T$especie[i]); m <- T$miniatura[!is.na(T$miniatura) & sub("_.*", "", T$especie) == g]
  if (length(m)) T$miniatura[i] <- m[1]
}
Q <- matrix(c(0, 1, 0, 0), 2, 2, dimnames = list(c("ausente", "presente"), c("ausente", "presente")))
analizar <- function(tr, x, mini, nsim = 50, pagel = TRUE) {
  H <- max(node.depth.edgelength(tr)); n_perd <- 0; clados <- c()
  for (n in c(1:Ntip(tr), (Ntip(tr) + 2):(Ntip(tr) + tr$Nnode))) {
    h <- if (n <= Ntip(tr)) tr$tip.label[n] else extract.clade(tr, n)$tip.label
    if (all(x[h] == "ausente")) {
      p <- tr$edge[tr$edge[, 2] == n, 1]
      if (!all(x[extract.clade(tr, p)$tip.label] == "ausente")) { n_perd <- n_perd + 1
        clados <- c(clados, paste(sort(unique(sub("_.*", "", h))), collapse = "+")) }
    }
  }
  sm <- make.simmap(tr, x, model = Q, nsim = nsim, pi = c(0, 1), message = FALSE)
  nl <- sapply(sm, function(s) sum(sapply(s$maps, function(m) length(m) > 1 && names(m)[1] == "presente")))
  pg <- c(NA, NA)
  if (pagel) for (k in 1:2) {
    y <- if (k == 1) ifelse(mini == "progenetic-miniature", "mini", "no") else
           ifelse(mini %in% c("progenetic-miniature", "proportioned-dwarf"), "mini", "no")
    ok <- !is.na(y); tt <- keep.tip(tr, tr$tip.label[ok])
    f <- tryCatch(fitPagel(tt, x[tt$tip.label], setNames(y[ok], tr$tip.label[ok])[tt$tip.label]), error = function(e) NULL)
    if (!is.null(f)) pg[k] <- f$P
  }
  list(dollo = n_perd, clados = paste(sort(clados), collapse = "; "), simmap = median(nl), simmap_rango = paste(range(nl), collapse = "-"),
       p_prog = pg[1], p_prog_enana = pg[2])
}
res <- list()
## (A) Fish Tree of Life
ft <- read.tree(file.path(D, "fishtree/full_20arboles.tre"))
sp <- unique(T$especie); sinon <- c(Danio_choprae = "Danio_choprai")
SIN_SP <- Sys.getenv("SIN_SP") == "1"   # variante: excluir especimenes no descritos (_sp), cuya posicion en FToL es solo por nombre
if (SIN_SP) sp <- sp[!grepl("_sp$", sp)]
SUF <- if (SIN_SP) "_sinSP" else ""
for (i in seq_along(ft)) {
  tr <- ft[[i]]; usados <- c(); mapa <- c(); sust <- c()
  for (s in sp) {   # 1.a pasada: nombres exactos (para que un congenere no ocupe la punta de otra especie del censo)
    s2 <- ifelse(s %in% names(sinon), sinon[s], s)
    if (s2 %in% tr$tip.label) { mapa[s] <- s2; usados <- c(usados, s2) }
  }
  for (s in setdiff(sp, names(mapa))) {   # 2.a pasada: congenere no usado
    cand <- sort(setdiff(grep(paste0("^", sub("_.*", "", s), "_"), tr$tip.label, value = TRUE), usados))
    if (length(cand)) { mapa[s] <- cand[1]; usados <- c(usados, cand[1]); sust <- c(sust, paste0(s, "->", cand[1])) }
  }
  tr <- keep.tip(tr, unname(mapa))
  est <- setNames(T$estado[match(names(mapa), T$especie)], mapa)[tr$tip.label]
  mini <- setNames(T$miniatura[match(names(mapa), T$especie)], mapa)[tr$tip.label]
  r <- analizar(tr, est, mini, nsim = 50, pagel = (i <= 5))
  res[[length(res) + 1]] <- data.frame(arbol = sprintf("FToL_completo_%02d", i), puntas = Ntip(tr), r)
  if (i == 1) { cat(sprintf("FToL: %d especies emparejadas, %d por congenere, %d sin genero en el arbol\n",
                            length(mapa), length(sust), length(sp) - length(mapa)))
                writeLines(sust, file.path(O, paste0("robustez_FToL_sustituciones", SUF, ".txt"))) }
  cat(sprintf("  arbol %02d: Dollo=%d simmap=%g (%s) p_prog=%s p_prog_enana=%s | %s\n", i, r$dollo, r$simmap, r$simmap_rango,
              format(r$p_prog, digits = 3), format(r$p_prog_enana, digits = 3), r$clados))
}
## (B) nuestro arbol de 24 genomas
tr <- read.tree(file.path(D, "arbol_ext/mf10.treefile"))
tr <- root(tr, outgroup = c("Astyanax_mexicanus", "Ictalurus_furcatus"), resolve.root = TRUE)
tr <- drop.tip(tr, c("Astyanax_mexicanus", "Ictalurus_furcatus"))
tr <- chronos(tr, lambda = 1, quiet = TRUE); class(tr) <- "phylo"
g <- T[T$fuente == "genoma", ]; est <- setNames(g$estado, g$especie)[tr$tip.label]
mini <- setNames(g$miniatura, g$especie)[tr$tip.label]
r <- analizar(tr, est, mini, nsim = 100)
res[[length(res) + 1]] <- data.frame(arbol = "ML_24_genomas", puntas = Ntip(tr), r)
cat(sprintf("  ML 24 genomas: Dollo=%d simmap=%g (%s) p_prog=%.3g p_prog_enana=%.3g | %s\n", r$dollo, r$simmap,
            r$simmap_rango, r$p_prog, r$p_prog_enana, r$clados))
R <- do.call(rbind, res); write.table(R, file.path(O, paste0("robustez_arboles", SUF, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nFToL (20 arboles): Dollo = %s perdidas\n", paste(names(table(R$dollo[grepl("FToL", R$arbol)])), collapse = "/")))
cat("FIN\n")
