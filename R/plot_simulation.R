#' @title Plots the SIR epidemic from the simulation
#'
#' @description Plots the total fraction of Susceptible, Infected and Recovered individuals at each time step.
#'
#' @param results The object containing the results from the simulation.
#' @param num_areas The number of areas in the simulation. Easy way to get it for your current simulation is
#'                  to use the dimension of your transmission kernel (e.g. if you have 4 age groups:
#'                  dim(expanded_D)[1]/4).
#' @param num_ages The number of ages in the simulation.
#' @param by_age Logical. If TRUE, will create a plot for each age category (only works with pretty = FALSE).
#' @param step Timestep used for the simulation, in days (default is 1 day).
#' @param pretty Logical. If TRUE, will render using ggplot (longer but nicer).
#' @param title Title to use for the plot.
#'
#' @return Creates a plot.
#'
#' @details The Infected category has a separate y-axis to make visualisation easier, due to the relatively low
#'          fraction of Infected at any given time compared to Susceptible and Recovered.
#'
#' @examples
#'
#' #Create a spatial dataset:
#' test_data = raster(nrow=10, ncol=10, xmn=1, xmx=100000, ymn=1, ymx=100000)
#' values(test_data) = runif(100, 1, 1000)
#'
#' #Calculate the parameters for the simulation:
#' prep_simulation(test_data)
#'
#' #Run the simulation:
#' results = run_simulation(test_data, expanded_D, contact_mat, beta)
#'
#' plot_simulation(results, dim(expanded_D)[1]/4, 4, pretty=T, title="SIR plot")
#'
#'
#' @export


plot_simulation = function(results, num_areas, num_ages, by_age=F, step=1, pretty=F, title=NULL){

  #inverse step for cleaner calculations involving this value below:
  step = 1/step

  #using ggplot:
  if(pretty == T){

    total = sum(results[1,2:dim(results)[2]])

    data = data.frame(Time = c(0:(max(results[,1])/step)),
                      S = rowSums(results[seq(1,max(results[,1])+1, step),2:(num_areas*num_ages+1)])/total,
                      I = rowSums(results[seq(1,max(results[,1])+1, step),(num_areas*num_ages+2):(num_areas*num_ages*2+1)])/total,
                      R = rowSums(results[seq(1,max(results[,1])+1, step),(num_areas*num_ages*2+2):dim(results)[2]])/total)

    ggplot()+
      geom_line(aes(x=data$Time, y=data$S, col="blue"))+
      geom_line(aes(x=data$Time, y=data$I*15, col="red"))+
      geom_line(aes(x=data$Time, y=data$R, col="green"))+
      scale_y_continuous(name = "Susceptible, Recovered", limits=c(0,1), breaks=seq(0,1,0.2), sec.axis = sec_axis(~./15 , name = "Infected", breaks=seq(0,1/15,0.01)))+
      scale_x_continuous(name = "Time (days)", breaks = seq(0,max(data$Time),20), minor_breaks = NULL)+
      scale_color_manual(name=NULL,breaks=c("blue","red","green"),values=c("royalblue","green3","red3"),labels=c("Susceptible","Infected", "Recovered"))+
      theme_bw()+
      theme(legend.box.background = element_rect(colour = "black"), legend.background = element_blank(),
            legend.justification=c(1,1), legend.position=c(1,1))+
      ggtitle(title)


    #default plotting:
  } else {


    if(by_age == FALSE){

      par(mar=c(5,4,4,5)+.1)

      total = sum(results[1,2:dim(results)[2]])

      plot(c(0:(max(results[,1])/step)), rowSums(results[seq(1,max(results[,1])+1, step),2:(num_areas*num_ages+1)])/total, col="blue", type="l", xlab="Time", ylab="S, R", ylim=c(0, 1), main=title)
      lines(c(0:(max(results[,1])/step)), rowSums(results[seq(1,max(results[,1])+1, step),(num_areas*num_ages*2+2):dim(results)[2]])/total, col="green")

      par(new=T)
      plot(c(0:(max(results[,1])/step)), rowSums(results[seq(1,max(results[,1])+1, step),(num_areas*num_ages+2):(num_areas*num_ages*2+1)])/total, type="l", col="red", xlab="",ylab="", xaxt="n",yaxt="n")
      axis(4)
      mtext("I",side=4,line=3)

      legend("topright", col=c("blue", "red", "green"), legend=c("S", "I", "R"), lty=1)

    } else {

      par(mfrow=c(2,2))

      S_1 = rowSums(results[,2:(num_areas+1)])
      S_2 = rowSums(results[,(num_areas+2):(num_areas*2+1)])
      S_3 = rowSums(results[,(num_areas*2+2):(num_areas*3+1)])
      S_4 = rowSums(results[,(num_areas*3+2):(num_areas*4+1)])

      I_1 = rowSums(results[,(num_areas*4+2):(num_areas*5+1)])
      I_2 = rowSums(results[,(num_areas*5+2):(num_areas*6+1)])
      I_3 = rowSums(results[,(num_areas*6+2):(num_areas*7+1)])
      I_4 = rowSums(results[,(num_areas*7+2):(num_areas*8+1)])

      R_1 = rowSums(results[,(num_areas*8+2):(num_areas*9+1)])
      R_2 = rowSums(results[,(num_areas*9+2):(num_areas*10+1)])
      R_3 = rowSums(results[,(num_areas*10+2):(num_areas*11+1)])
      R_4 = rowSums(results[,(num_areas*11+2):(num_areas*12+1)])

      total1 = sum(S_1[1],I_1[1],R_1[1])
      total2 = sum(S_2[1],I_2[1],R_2[1])
      total3 = sum(S_3[1],I_3[1],R_3[1])
      total4 = sum(S_4[1],I_4[1],R_4[1])


      plot(results[,1], S_1/total1, type="l", col="blue", ylab="S,I,R fraction", xlab="Time", main="0-4 years old", ylim=c(0, 1))
      lines(results[,1], I_1/total1, col="red")
      lines(results[,1], R_1/total1, col="green")
      legend("topright", col=c("blue", "red", "green"), legend=c("S", "I", "R"), lty=1)


      plot(results[,1], S_2/total2, type="l", col="blue", ylab="S,I,R fraction", xlab="Time", main="5-19 years old", ylim=c(0, 1))
      lines(results[,1], I_2/total2, col="red")
      lines(results[,1], R_2/total2, col="green")
      legend("topright", col=c("blue", "red", "green"), legend=c("S", "I", "R"), lty=1)


      plot(results[,1], S_3/total3, type="l", col="blue", ylab="S,I,R fraction", xlab="Time", main="20-64 years old", ylim=c(0, 1))
      lines(results[,1], I_3/total3, col="red")
      lines(results[,1], R_3/total3, col="green")
      legend("topright", col=c("blue", "red", "green"), legend=c("S", "I", "R"), lty=1)


      plot(results[,1], S_4/total4, type="l", col="blue", ylab="S,I,R fraction", xlab="Time", main="65+ years old", ylim=c(0, 1))
      lines(results[,1], I_4/total4, col="red")
      lines(results[,1], R_4/total4, col="green")
      legend("topright", col=c("blue", "red", "green"), legend=c("S", "I", "R"), lty=1)


    }

    #reset graphics parameters
    par(mfrow=c(1,1))

  }

}

