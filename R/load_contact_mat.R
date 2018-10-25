#' @title Creates the age mixing matrix
#'
#' @description Creates the age mixing matrix for four age groups: 0-4, 5-19, 20-64 and 65+, or for no age groups.
#'
#' @param age Logical. If FALSE, will assume homogeneous mixing between age groups in the population and therefore
#'            return a matrix to enable this in the other functions (default is TRUE).
#'
#' @return Returns one matrix object.
#'
#' @examples
#'
#' contact_mat = load_contact_mat()
#' View(contact_mat)
#'
#' @export

load_contact_mat = function(age=T){

  if(age == TRUE){

    contact_mat = matrix(c(37.4622640266199,13.2337799407673,9.35866526693108,5.27807222067513,
                           17.2304141889828,98.1983003738366,17.0186152145963,10.1131975048866,
                           9.46784315245156,9.4416088929148,16.22285757548,5.7675253611147,
                           1.38284918679668,1.26680284573205,1.08367881504336,3.88324564380799),
                         ncol=4, nrow=4, byrow=T,
                         dimnames = list(c("0-4","5-19","20-64","65+"),
                                         c("0-4","5-19","20-64","65+")))

  } else {

    contact_mat = matrix(1, 1, 1)

  }

  return(contact_mat)

}
