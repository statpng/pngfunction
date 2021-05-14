https://r-graphics.org/recipe-output-fonts-pdf

library(extrafont)

p <- ggplot2::ggplot(iris, aes(Sepal.Length, Sepal.Width, color=Species)) + geom_point()
ggsave(p, filename = "tmp.pdf")

embed_fonts("tmp.pdf")
