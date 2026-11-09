# Paleta derivada de Leon et al. 2024, Mol Biol Evol 41(9):msae166 (CC BY 4.0); HEX medidos en PMC11371425.
leon <- c(presente = "#664BBC", ausente = "#FA41CA", parcial = "#FEDC00", externo = "#71D226",
          sintenia = "#66FEFF", hueco = "#E5C0DA",
          ctrl1 = "#EBEBEB", ctrl2 = "#C4BAB9", ctrl3 = "#000000",
          texto = "#000000", texto2 = "#525252", rejilla = "#EBEBEB")
leon_secuencial <- c("#FDE2CD", "#F0ACA9", "#FE0000", "#9A0000")   # solo para variables continuas
library(ggplot2)
theme_leon <- function(base_size = 8, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(text = element_text(colour = leon[["texto"]]),
          axis.text = element_text(colour = leon[["texto"]]),
          axis.line = element_line(colour = leon[["texto2"]], linewidth = 0.3),
          panel.grid.major = element_line(colour = leon[["rejilla"]], linewidth = 0.3),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.key = element_blank())
}
# Estado de DCX: color Y forma (nunca solo color)
scale_colour_dcx <- function(...) scale_colour_manual(values = c(presente = leon[["presente"]], ausente = leon[["ausente"]]), ...)
scale_shape_dcx  <- function(...) scale_shape_manual(values = c(presente = 16, ausente = 21), ...)   # lleno / vacio
scale_fill_leon_c <- function(...) scale_fill_gradientn(colours = leon_secuencial, ...)
# Ejemplo:
# ggplot(df, aes(x, y, colour = estado, shape = estado)) + geom_point(size = 2, stroke = 1, fill = "white") +
#   scale_colour_dcx() + scale_shape_dcx() + theme_leon()
# ggsave("fig.pdf", width = 120, height = 100, units = "mm", device = cairo_pdf)
# ggsave("fig.png", width = 120, height = 100, units = "mm", dpi = 600)
