#Create likelihood function 
gamma_likelihood = function(pars,se,pi_ref,mbar = 1) {
  #Number of datapoints, is length of o_sd
  n      = length(se)
  #Set variables
  c      = 1 / pars[1] - 1
  alpha  = c * ( pi_ref )
  beta   = ( 1 - pi_ref ) * c
  
  #Define the marginal likelihood function (which will we then integrate)
  marginal = function(z, 
                      se,
                      alpha, 
                      beta,
                      noise){
    log_result = log( 1 / ( ( z^-1*se*mbar - 4 )^.5 * z^-1*se*mbar ^ ( 3 / 2 ) ) ) +
      log( ( dbeta( ( 0.5 - ( ( z^-1*se*mbar - 4 ) ^ .5 ) / ( 2 * ( z^-1*se*mbar ) ^ .5 ) ) , alpha, beta ) +  
               dbeta( ( 0.5 + ( ( z^-1*se*mbar - 4 ) ^ .5 ) / ( 2 * ( z^-1*se*mbar ) ^ .5 ) ) , alpha, beta ) ) ) + 
      dgamma( z, shape = 1 / pars[2], scale = pars[2], log = TRUE) +
      log( z^-1 )
    
    return(exp(log_result))
  }
  
  Likelihood = c()
  #Calculate the likelihood across each observation
  for (i in 1:n) {
    #Integrate out this nuissance term
    #Create a check to figure out where the lower bound should be
    lower_bound = 0
    while (length(Likelihood) == (i-1)) {
      dat = try(expr = integrate(f = marginal, lower = lower_bound, upper = se[i]*mbar/4, se = se[i], alpha = alpha[i], beta = beta[i]),silent = TRUE)
      if(class(dat) != "try-error"){
        Likelihood[i] = integrate(f = marginal, lower = lower_bound, upper = se[i]*mbar/4, se = se[i], alpha = alpha[i], beta = beta[i])$value
      }
      else{
        lower_bound = lower_bound + 0.01}
    }
  }
  #print(Likelihood)
  Likelihood = Likelihood[!is.na(Likelihood) & !is.infinite(Likelihood)]
  return(-sum(log(Likelihood)))
}