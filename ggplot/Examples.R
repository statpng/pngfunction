# In a histogram, highlight specific samples from the total samples.

library(ggplot2)
dat <-rnorm(80)
dat <-data.frame(dat)
p <- ggplot(dat, aes(x=dat))+geom_histogram()
p
## filtering...특정 데이터만 추출하기
dat %>%
  filter(dat >= -.5 & dat <= .5) -> dat_filtered

p + geom_rug(data = dat_filtered, aes(x = dat), colour="blue",  inherit.aes = F)



# Hangul Font
# https://r-graphics.org/recipe-output-fonts-pdf

library(extrafont)
p <- ggplot2::ggplot(iris, aes(Sepal.Length, Sepal.Width, color=Species)) + geom_point() ggsave(p, filename = "tmp.pdf")
embed_fonts("tmp.pdf")


