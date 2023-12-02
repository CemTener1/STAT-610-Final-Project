#Import the dataset
H3N2_78 <- matrix(c(66, 87, 25, 22, 4, 13, 14,15, 9, 4,NA, 4,4,9,1,NA,NA, 4,3,1,
                 NA,NA,NA,1,1,NA,NA,NA,NA,0), ncol = 5, nrow = 6, byrow =TRUE)

H3N2_81 <- matrix(c(44, 62, 47, 38, 9, 10,13, 8, 11,5, NA,9,2,7,3, NA,NA, 3,5,1,
                    NA,NA,NA, 1,0,NA,NA,NA,NA,1 ), ncol = 5, nrow = 6, byrow =TRUE)


H1N1_76 <- matrix(c(9, 12, 18, 9, 4, 1,6,6,4,3,NA,2,3,4,0,NA,NA, 1,3,2,NA,NA,NA,0,
                    0, NA,NA,NA,NA,0), nrow = 6,ncol = 5, byrow = TRUE)

H1N1_79 <-matrix(c(15,12,4, 11,17,4,NA,21,4,NA,NA, 5), nrow = 4, ncol = 3, byrow = TRUE)


#Calculate the distance using Frobenious norm
## ||D1 - D*|| Frobenious norm is the square root of the sum of the squares of element-wise difference 
## ||D1 - D*||F+ negative elements (differences) are replaced with 0

distance <- function(D1,D2, D_star1, D_star2){
  return(0.5*(norm(pmax(D1 - D_star1), type = "F")*norm(D2 - D_star2, type = "F"))) 
  
}


#Initialize the prior and 
N = 10000
qc1 = runif(N)
qh1 = runif(N)
qc2 = runif(N)
qh2 = runif(N)
  
  
generate_abc <- function(observed, qc, qh, epsilon,N, househould_size, epsilon){
 accepted_sample = data.frame(qc1 = numeric(0),
                              qh1 = numeric(0),
                              qc2 = numeric(0),
                              qh2 = numeric(0))
 
 for (i in 1:N){
   #simulate data for given parameters
   simulated_data <- simulate(qc1[i], qh1[i], qc2[i], qh2[i], househould_size)
   #assign the simulated data to a new variable 
   D_star1 <- simulated_data$D_star1
   D_star2 <- simulated_data$D_star2
   #calculate the distance
   dist <- distance(D1, D2, D_star1, D_star2)
   
   #Accept/Reject Condition
   if (dist < epsilon){
    accepted_sample[nrow(accepted_sample) + 1,] <-list(qc1, qh1, qc2, qh2)
   }
 }
  
}




































































