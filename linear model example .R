covariate <- sample(0:1, 100, replace=TRUE)
exposure  <- runif(100,0,1)+(0.3*covariate)
outcome   <- 2.0+(0.5*exposure)+(0.25*covariate)