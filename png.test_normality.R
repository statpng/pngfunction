png.test_normality <- function(x, test=c("shapiro", "ks"), bestNormalize=FALSE){

  name.x <- deparse(substitute(x))
  test <- match.arg(test)
  
  if( bestNormalize ){
    library(bestNormalize)

    func_list <- c("arcsinh_x(x)",
                   "boxcox(x)",
                   "log_x(x)",
                   "x",
                   "orderNorm(x)",
                   "sqrt_x(x)",
                   "yeojohnson(x)")
    tranformed_x <- sapply( func_list, function(y) {
      res <- eval(parse(text=y))
      if( is.vector(eval(parse(text=y))) ){
        res
      } else {
        res$x.t
      }
      })
  } else {
  func_list <- c("1/x^2",
                 "1/x",
                 "1/sqrt(x)",
                 "log(x)",
                 "sqrt(x)",
                 "x",
                 "x^2",
                 "x^(-0.014)")
  tranformed_x <- sapply( func_list, function(y) eval(parse(text=y)) )
  }

  
  pvalue <- NULL
  for( j in 1:ncol(tranformed_x) ){
    
    xj <- tranformed_x[,j,drop=T]
    if( test == "shapiro" ){
      pvaluej <- shapiro.test(xj)$p.value
    } else if ( test == "ks" ){
      pvaluej <- ks.test(xj, "pnorm")$p.value
    }
    
    pvalue[j] <- max( pvaluej, 1e-16 )
  
  }
  
  library(ggplot2)
  p <- tidyr::gather(as.data.frame(tranformed_x), variable, value) %>% 
    mutate(variable=factor(variable, levels=func_list),
           pvalue=factor(variable, levels=func_list, labels=ifelse( pvalue<1e-3, "p < 0.001", paste0("p = ", round(pvalue,3)) ))) %>% 
    group_by(variable) %>% 
    mutate( density_y = dnorm( x=value, mean=mean(value, na.rm=TRUE), sd=sd(value, na.rm=TRUE))
            # range_x = max(value, na.rm = TRUE),
            # range_y = max(density_y, na.rm = TRUE)
            ) %>% 
    filter(!is.na(value)) %>% 
      ggplot(aes(x=value)) +
      geom_histogram(aes(y=..density..), fill = "lightyellow", color = "black") +
      geom_line(aes(value, density_y), color="red") +
      geom_text(aes( label=pvalue, x=Inf, y=Inf ), hjust="inward", vjust="inward" ) +
      facet_wrap(~variable, scales = "free") +
      theme_bw() +
    scale_x_continuous(name=NULL)+
      labs(title=paste0("Transformation of ", name.x), 
           subtitle = paste0("p values by ", test) )
  p
}

# hist_ordinary <- transformDensity(moonBook::acs$TG, test="shapiro", bestNormalize = FALSE)
# hist_best <- transformDensity(moonBook::acs$TG, test="shapiro", bestNormalize = TRUE)
# 
# hist_ordinary %>% ggsave(filename="test_ordinary.jpeg", width=10, height=5)
# hist_best %>% ggsave(filename="test_best.jpeg", width=10, height=5)
# 
# gridExtra::grid.arrange(
#   hist_ordinary, hist_best, nrow=1
# ) %>% ggsave(filename="test_combined.jpeg", width=20, height=5)
