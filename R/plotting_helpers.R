theme_mnp <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

save_plot <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(filename)
}
