

findcritical <- function() {
  #Find the critical value z such that P(Z>=z)=alpha 
  value <- 0 
  alpha <- readline("What is alpha over 2 (1-alpha is the confidence of the interval)")
  alpha <- as.numeric(alpha)
  print("Normal distribution, press 1")
  print("t distribution, press 2")
  print("Chi-squared distribution, press 3")
  distribution <- readline("What distribution do you want to use")
  distribution <- as.numeric(distribution)
  if (distribution == 1) {
    print("Standard Normal")
    value <- qnorm(1-alpha,mean=0,sd=1)
  } else if (distribution == 2) {
    df <- readline("What is the degree of freedom for the t distribution")
    df <- as.numeric(df)
    value <- qt(alpha, df, lower.tail = FALSE, log.p = FALSE)
    cat("t distribution with degree of freedom", df)
  } else if (distribution == 3) {
    df <- readline("What is the degree of freedom for the chi-squared distribution")
    df <- as.numeric(df)
    value <- qchisq(alpha, df, ncp = 0, lower.tail = FALSE, log.p = FALSE)
    cat("Chi-squared distribution with degree of freedom", df)
  }
  print("Warning. The z given is such that P(Z>=z)=alpha")
  return(value) 
}