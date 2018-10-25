#' @title Estimates beta from a given R0
#'
#' @description Estimates beta for an epidemic with a given R0 in a specific population in a spatial dataset,
#'              and calculates the expanded kernel matrix (expanded_D) to use when simulating an epidemic.
#'
#' @return Returns the value of beta and assigns the expanded_D matrix to the global environment.
#'
#' @param spatial_data The spatial_dataset containing the population data.
#' @param dist_kernel The distance kernel matrix.
#' @param contact_mat The contact matrix for mixing between age groups.
#' @param R0 The desired value for R0.
#' @param sigma The desired value for the recovery rate.
#'
#' @details This function is automatically executed when using the \code{\link{prep_simulation}} function. It uses
#'          the Next Generation Matrix approach to derive beta from R0.
#'
#' @examples
#'
#' #Create a spatial dataset:
#' test_data = raster(nrow=10, ncol=10, xmn=1, xmx=100000, ymn=1, ymx=100000)
#' values(test_data) = runif(100, 1, 1000)
#'
#' #Calculate distance kernel matrix and load age mixing matrix:
#' dist_mat = calc_dist_mat(test_data)
#' dist_kernel = calc_dist_kernel(dist_mat, dist_c = 87, test_data, alpha=0.95, p=6.6, p2=1.53, aa=35)
#' load_contact_mat()
#'
#' beta = calc_beta(test_data, dist_kernel, contact_mat, R0=1.8, sigma=1/2.6)
#'
#' @export


calc_beta = function(spatial_data, dist_kernel, contact_mat, R0=1.8, sigma=1/2.6){

  if(class(spatial_data) == "RasterLayer"){

    populated_areas = which(!is.na(spatial_data@data@values))

    N = spatial_data@data@values[populated_areas]

  } else if(class(spatial_data) == "SpatialPolygonsDataFrame"){

    N = spatial_data@data$population

  } else{

    stop("Incorrect spatial dataset! The dataset must be either a RasterLayer or a SpatialPolygonsDataFrame.")

  }


  num_areas = dim(dist_kernel)[1]   #derive number of areas from distance kernel
  num_ages = dim(contact_mat)[1]    #derive number of age categories from contact matrix


  #!!! work in progress, only supports 4 age categories or no age categories at all right now !!!#

  if(num_ages == 4){

    N = matrix(N, nrow=num_areas, ncol=num_ages)

    #set minimum population in an area to 1:
    N[which(N<1)] = 1

    NN0 = as.vector(N)

    N[,1] = N[,1]*(5/81)
    N[,2] = N[,2]*(14/81)
    N[,3] = N[,3]*(46/81)
    N[,4] = N[,4]*(16/81)
    #this way, i in the matrix is the area and j the age group e.g. N[1,1] gives pop 0-4 in area 1


  } else if(num_ages == 1){

    N[which(N<1)] = 1

    NN0 = as.vector(N)

  } else {

    stop("Unsupported number of age categories, currently only supports 4 or none.")

  }


  Sstart = matrix(N, ncol=1)
  Sstart = matrix(Sstart, ncol=nrow(Sstart), nrow=nrow(Sstart))

  K1 = kronecker(diag(num_ages), dist_kernel)

  NNbar = matrix(N, ncol=1)
  Kbar = kronecker(matrix(1,num_ages,num_ages), dist_kernel)

  Mj=t(Kbar)%*%NNbar
  Mj[which(Mj==0)]=1
  Mjover=1/Mj
  Mjover = t(Mjover)
  Mjover = matrix(rep(Mjover, num_areas*num_ages), nrow=num_areas*num_ages, byrow=T)

  keye = diag(num_areas)
  kxeye = matrix(1, nrow=num_areas, ncol=num_areas) - keye
  Cbar = kronecker(contact_mat, keye) + kronecker(contact_mat, kxeye)

  DD=(Sstart*K1*Mjover)%*%(t(Kbar)*Cbar)

  X=DD/sigma


  ### this extra bit just calculates D, necessary for FOI calculation, it's then faster to run the simulation
  Ni=matrix(NN0,nrow=length(NN0),ncol=length(NN0),byrow=F)
  Nj=t(Ni)

  proper_D = (K1*Mjover)%*%(t(Kbar)*Cbar)
  proper_D = proper_D*Nj

  assign("expanded_D", proper_D, envir=.GlobalEnv)

  ###


  eigen_vals = eigen(X, symmetric=F, only.values = T)$values

  R0a = max(Re(eigen_vals))

  beta=R0/R0a

  return(beta)

}
