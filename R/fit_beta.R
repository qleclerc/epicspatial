#' @title Fits the incidence from one resolution to another via beta
#'
#' @description Adjusts beta to fit the incidence of one resolution onto another.
#'
#' @param fit_data The incidence data to use for fitting.
#' @param spatial_data The spatial dataset to be fitted.
#' @param interval The interval on which the fit should be performed, in days (default is c(1,100) - i.e. the fit
#'                 will be performed on the first 100 days of the epidemic).
#' @param stoch Logical. If TRUE (default), fitting will be performed using the median of x stochastic runs. If
#'              FALSE, fitting will be performed using the deterministic output.
#' @param num_runs Number of stochastic runs to extract the median incidence for one beta value (default is 100).
#'                 Higher values will be longer to run, but yield more accurate results. (only considered if
#'                 stoch = TRUE)
#' @param search_range Search range for the fitted value of beta (default is c(0.02,0.05)).
#'
#' @details The fitting is performed by assuming the fitted incidence belongs to a Poisson distribution with the
#'          expected values being the values from fit_data. Consequently, the values from fit_data must be integers.
#'
#' @return Returns the fitted value of beta.
#'
#' @examples
#'
#' #Create a high resolution spatial dataset:
#' htest_data = raster(nrow=20, ncol=20, xmn=1, xmx=100000, ymn=1, ymx=100000)
#' values(htest_data) = runif(400, 1, 1000)
#'
#' #Calculate the median incidence for the high resolution data:
#' prep_simulation(htest_data)
#' results_high = run_multi_stoch(100, htest_data, expanded_D, contact_mat, beta, t_max=100)
#' results_high = round(results_high$Median)
#'
#'
#' #Create a low resolution spatial dataset (4 times less areas than htest_data):
#' ltest_data = aggregate(htest_data, fact=2, fun=sum)
#'
#' #Fit low resolution on high resolution:
#' fitted_beta = fit_beta(results_high, ltest_data, interval=c(1,100))
#'
#' @export


fit_beta = function(fit_data, spatial_data, interval=c(1,100), stoch=TRUE, num_runs=100, search_range=c(0.02,0.05)){

  if(length(fit_data) != (interval[2]-interval[1]+2)) stop("Fit data length does not match the fitting interval!")

  prep_simulation(spatial_data)

  #log likelihood calculation:
  find_beta = function(fit_data, test_beta, stoch, spatial_data){

    if(stoch==TRUE){
      results = run_multi_stoch(num_runs, spatial_data, expanded_D, contact_mat, test_beta, t_max=interval[2])
      results = results$Median[interval[1]:length(results$Median)]
    }

    if(stoch==FALSE){
      results = run_simulation(spatial_data, expanded_D, contact_mat, test_beta, t_max=interval[2])
      num_ages = 4
      num_areas = dim(expanded_D)[1]/num_ages
      results = rowSums(results[,2:(num_areas*num_ages+1)])
      results = c(0,abs(diff(results)))
      results = results[interval[1]:length(results)]
    }

    #plot the fitted curve against the target each time a beta is estimated to visualise progress:
    plot(x=seq(interval[1],interval[2]+1,1), y=fit_data, type="l", main=paste0("Beta: ", test_beta))
    lines(x=seq(interval[1],interval[2]+1,1), y=results, col="red")
    legend("topright", legend=c("Target", "Fitted"), col=c("black","red"), lty=1)


    #calculate total log likelihood of fit:
    ll = sum(dpois(fit_data, results, log=T))

    return(ll)

  }

  #optimise function to find beta yielding highest total log likelihood:
  best_beta = optimise(function(test_beta) find_beta(fit_data, test_beta, stoch, spatial_data), search_range, maximum = T)

  return(best_beta$maximum)

}
