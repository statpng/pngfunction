library(ggplot2)
dat <-rnorm(80)
dat <-data.frame(dat)
p <- ggplot(dat, aes(x=dat))+geom_histogram()
p
## filtering...특정 데이터만 추출하기
dat %>%
  filter(dat >= -.5 & dat <= .5) -> dat_filtered

p + geom_rug(data = dat_filtered, aes(x = dat), colour="blue",  inherit.aes = F)
