#' @title Performs pre-calculations necessary for the simulation of an epidemic
#'
#' @description Creates the age mixing matrix and the spatial kernel using the \code{\link{load_contact_mat}},
#'              \code{\link{calc_dist_mat}}, \code{\link{calc_dist_kernel}} and \code{\link{calc_beta}} functions.
#'
#' @param spatial_data The spatial dataset containing the population data.
#' @param R0 The desired R0 for the epidemic
#' @param sigma The desired recovery rate for the epidemic
#' @param age Logical. If TRUE, the simulation will assume four distincts age groups in the population (default).
#'            If FALSE, the simulation will assume homogeneous mixing between age groups in the population.
#' @param k_params A list containing the parameter values to be used for the kernel calculation. If omitted,
#'                 will use the following default values: dist_c=87; alpha=0.95; p=6.6; p2=1.53; aa=35
#'
#' @return Returns two matrix objects and one value for beta.
#'
#' @examples
#'
#' #Create a spatial dataset:
#' test_data = raster(nrow=10, ncol=10, xmn=1, xmx=100000, ymn=1, ymx=100000)
#' values(test_data) = runif(100, 1, 1000)
#'
#' prep_simulation(test_data)
#'
#'
#' @export


prep_simulation = function(spatial_data, R0 = 1.8, sigma= 1/2.6, age=TRUE, k_params=NA){

  assign("contact_mat", load_contact_mat(age), envir=.GlobalEnv)

  print("Calculating distance matrix...")

  d = calc_dist_mat(spatial_data)

  print("Calculating kernel matrix...")

  if(is.na(k_params)) k_params = c("dist_c"=87, "alpha"=0.95, "p"=6.6, "p2"=1.53, "aa"=35)

  #Must change this so that kernel parameters can be modified easily outside function!
  kernel = calc_dist_kernel(d, dist_c=k_params["dist_c"], spatial_data, alpha=k_params["alpha"], p=k_params["p"],
                            p2=k_params["p2"], aa=k_params["aa"])

  print("Calculating beta...")

  assign("beta", calc_beta(spatial_data, kernel, contact_mat, R0=R0, sigma=sigma), envir=.GlobalEnv)


}
