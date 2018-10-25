#' @title Calculates the spatial kernel values
#'
#' @description Calculates the spatial kernel values between all areas in a spatial dataset using their distance and
#'              population. This uses a gravity kernel structure (either matched or smooth depending on the threshold
#'              distance).
#'
#' @param dist_mat The distance matrix between all areas.
#' @param dist_c The threshold distance for the matched kernel. Note that setting this to an extremely high value
#'               will effectively make the gravity kernel a smooth kernel with a single function since the threshold
#'               distance will never be crossed.
#' @param spatial_data The spatial dataset containing the population data
#' @param alpha The destination population power.
#' @param p The distance power.
#' @param p2 The secondary distance power for the matched kernel.
#' @param aa The offset distance.
#' @param delta The home cell modifier (to compensate the distance within cells equal to 0).
#'
#' @return Returns one matrix object containing the spatial kernel values.
#'
#' @examples
#'
#' #Create a spatial dataset:
#' test_data = raster(nrow=10, ncol=10, xmn=1, xmx=100000, ymn=1, ymx=100000)
#' values(test_data) = runif(100, 1, 1000)
#'
#' #Calculate distance matrix
#' dist_mat = calc_dist_mat(test_data)
#'
#' dist_kernel = calc_dist_kernel(dist_mat, dist_c = 87, test_data, alpha=0.95, p=6.6, p2=1.53, aa=35)
#'
#' @export


calc_dist_kernel = function(dist_mat, dist_c, spatial_data, alpha, p, p2, aa, delta=0.3){

  if(class(spatial_data) == "RasterLayer"){

    populated_areas = which(!is.na(spatial_data@data@values))

    N = spatial_data@data@values[populated_areas]

  } else if(class(spatial_data) == "SpatialPolygonsDataFrame"){

    N = spatial_data@data$population

  } else{

    stop("Incorrect spatial dataset! The dataset must be either a RasterLayer or a SpatialPolygonsDataFrame.")

  }


  N[which(N<1)] = 1

  #set up matrix to fill in:
  dist_kernel = matrix(0, nrow=nrow(dist_mat), ncol=ncol(dist_mat))

  #value at threshold distance to ensure a continuous function:
  A = (dist_c^p2)/((1+dist_c/aa)^(p))

  #pre-threshold values:
  K = (1+dist_mat/aa)^(p)

  #only adds delta if i == j (i.e. if calculating transmission kernel within an area)
  delta_diag = matrix(1, nrow=nrow(dist_mat), ncol=ncol(dist_mat))
  diag(delta_diag) = 1+delta
  K = K*delta_diag

  #post-threshold values:
  K2 = A/(dist_mat^p2)

  #CHANGE THIS: useless calculations being performed, currently calculating K and K2 for ALL distances, better to
  #only calculate the relevant one based on the distance value


  for (i in 1:length(N)) {

    for (j in 1:length(N)) {

      if(dist_mat[i,j] > dist_c){

        dist_kernel[i,j] = ((N[j]^alpha))*(K2[i,j])

      } else {

        dist_kernel[i,j] = ((N[j]^alpha))/(K[i,j])

      }

    }

  }


  #normalise so that rowSums = 1:
  dist_kernel = dist_kernel/rowSums(dist_kernel)

  return(dist_kernel)

}

