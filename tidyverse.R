capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}
global_labeller <- labeller(
  vore = capitalize,
  # conservation = conservation_status,
  conservation2 = label_wrap_gen(10),
  .default = label_both
)
# usage: facet_grid( Group1~Group2, labeller = global_labeller )


# For library(patchwork)
label_plot <- function(label, angle) {
    # usage: label_plot("LA", 0)
    ggplot() + geom_text(aes(x = 0, y = 0, label = label), size = 6, fontface = 2, angle=angle) + theme_void()
  }
  
  
