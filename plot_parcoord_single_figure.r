library(MASS)
# /////////////////////////////////////////

generate_parcoord_plot_with_costs<- function(points_dataframe, fraction_of_maxforce_for_task, my_column_names, fmax_vector){
  points_with_costs <- fill_costs(points_dataframe, 7, fmax_vector)
  points_with_costs$TaskForce <- rep(fraction_of_maxforce_for_task, length(points_with_costs[,1]))
  colnames(points_with_costs) <- my_column_names
  
  # parcoord(points_with_costs, var.label=TRUE, ylim=c(0,1))
  
  nonweighted_max_observed_cost <- max(points_with_costs$L1)
  weighted_max_observed_cost <- max(points_with_costs$L1W)
  
  #unweighted
  points_with_costs$L1 <- points_with_costs$L1 / nonweighted_max_observed_cost
  points_with_costs$L2 <- points_with_costs$L2 / nonweighted_max_observed_cost
  #weighted
  points_with_costs$L1W <- points_with_costs$L1W / weighted_max_observed_cost
  points_with_costs$L2W <- points_with_costs$L2W / weighted_max_observed_cost
  
  require(GGally)
  require(ggplot2)
  p_raw <- ggparcoord(points_with_costs, scale='globalminmax', alpha=0.025, boxplot=FALSE, mapping=ggplot2::aes(colour="midnightblue"))
  p <- p_raw + theme_bw() + theme(
    panel.grid.major.x = element_line(color = "black", size = 0.5),
    panel.grid.major = element_blank(),
    legend.position = "none"
) + geom_hline(yintercept=seq(0, 1, by=0.2))
  return(p)
}

# fills raw cost columns (14-19 for finger)
# section 3 of dataframe
fill_costs = function(df,m, fmax_vector) {
  F0 = matrix(fmax_vector,dim(df)[1],m,byrow=TRUE)
  # L1
  df[,(m+1)] = rowSums(df[,1:m])
  # L2
  df[,(m+2)] = (rowSums(df[,1:m]^2))^(1/2)
  # Lw1
  df[,(m+3)] = rowSums(df[,1:m]*F0)
  # Lw2
  df[,(m+4)] = (rowSums((df[,1:m]*F0)^2))^(1/2)

  
  return(df)

}

frac_of_indexfinger_max <- function(task_force_newtons) {
  fraction_of_maxforce_for_task <- task_force_newtons/28.81155463796023
  return(fraction_of_maxforce_for_task)
}

my_column_names <- c("FDP",
                    "FDS",
                    "EIP",
                    "EDC",
                    "DI",
                    "PI",
                    "LUM",
                    "L1",
                    "L2",
                    "L1W",
                    "L2W",
                    "TaskForce")

fmax_vector=c(123, 219, 124.8, 129.6, 23.52, 21.6, 91.74) #newtons each muscle is maximally capable of producing.
#tasks chosen: c(0,.2,.4,.6,.8,.1)*28.811554635079077
points_a <- read.csv("n_1000_alphalen_1000/finger_forcevector_0.0_1484767835427.csv", header=FALSE)
points_b <- read.csv("n_1000_alphalen_1000/finger_forcevector_5.768079006021837_1484784224782.csv", header=FALSE)
points_c <- read.csv("n_1000_alphalen_1000/finger_forcevector_11.507317617013564_1484878583772.csv", header=FALSE)
points_d <- read.csv("n_1000_alphalen_1000/finger_forcevector_17.275396623035405_1484881413230.csv", header=FALSE)
points_e <- read.csv("n_1000_alphalen_1000/finger_forcevector_23.04347562905724_1484876827322.csv", header=FALSE)
points_f <- read.csv("n_1000_alphalen_1000/finger_forcevector_28.811554635079077_1484891275183.csv", header=FALSE)

plot_a <- generate_parcoord_plot_with_costs(points_a, frac_of_indexfinger_max(0.0), my_column_names, fmax_vector)
plot_b <- generate_parcoord_plot_with_costs(points_b, frac_of_indexfinger_max(5.768079006021837), my_column_names, fmax_vector)
plot_c <- generate_parcoord_plot_with_costs(points_c, frac_of_indexfinger_max(11.507317617013564), my_column_names, fmax_vector)
plot_d <- generate_parcoord_plot_with_costs(points_d, frac_of_indexfinger_max(17.275396623035405), my_column_names, fmax_vector)
plot_e <- generate_parcoord_plot_with_costs(points_e, frac_of_indexfinger_max(23.04347562905724), my_column_names, fmax_vector)
plot_f <- generate_parcoord_plot_with_costs(points_f, frac_of_indexfinger_max(28.811554635079077), my_column_names, fmax_vector)

require(gridExtra)
  combined_figure <- grid.arrange(
                  plot_a,
                  plot_b,
                  plot_c,
                  plot_d,
                  plot_e,
                  plot_f,
                  ncol=1)

ggsave("figure_5_parcoord.pdf", combined_figure, width=8, height=10)