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

long_df <- pivot_longer(accepted_samples, 
                        cols = c(qc1, qh1, qc2, qh2), 
                        names_to = c(".value", "group"), 
                        names_pattern = "(..)(.)")

# Plotting
ggplot(long_df, aes(x = qc, y = qh, color = group)) +
  geom_point() +
  labs(title = "Scatter Plot by Group",
       x = "QC Axis",
       y = "QH Axis") +
  theme_minimal()



########## Ky
household_size <- list(H3N2_78 <- colSums(H3N2_78,na.rm = T),
                       H3N2_81 <-colSums(H3N2_81,na.rm = T),
                       H1N1_76 <- colSums(H1N1_76,na.rm = T),
                       H1N1_79 <- colSums(H1N1_79,na.rm = T))
names(household_size) <- c('H3N2_78','H3N2_81','H1N1_76','H1N1_79')
calculate_w_matrix <- function(qc, qh, n) {
  # Initialize the matrix with NA values
  w_matrix <- matrix(NA, nrow = n+1, ncol = n)
  
  # Set the first row (w_0s)
  for (s in 1:n){
    w_matrix[1,s] <- qc^(s-1)
  }
  
  # Iterate over the rows
  for (j in 2:(n+1)) {
    # Calculate the diagonal element
    w_matrix[j, j-1] <- 1 - sum(w_matrix[1:(j-1), j-1], na.rm = T)
    
    # Fill the rest of the row
    for (s in (j-1):n) {
      if (s > j-1) {
        w_matrix[j, s] <- choose(s, j-1) * w_matrix[j, j-1] * (qc * qh^(j-1))^(s-j+1)
      }
    }
  }
  rownames(w_matrix) <- 0:n
  colnames(w_matrix) <- 1:n
  
  # Return the completed matrix
  return(w_matrix)
}

make_table_2 <- function(qc1,qh1,qc2,qh2){
  ## for table 2, both year period haves col num 5
  n1 <- 5
  n2 <- 5
  p1 <- calculate_w_matrix(qc1,qh1,n1)
  p2 <- calculate_w_matrix(qc2,qh2,n2)
  m1 <- p1 * household_size$H3N2_78
  m2 <- p2 * household_size$H3N2_81
  return(cbind(m1,m2))
}

make_table_3 <- function(qc1,qh1,qc2,qh2){
  ## for table 3, 75-76 has col 5 and 78-79 has col 3
  n1 <- 5
  n2 <- 3
  p1 <- calculate_w_matrix(qc1,qh1,n1)
  p2 <- calculate_w_matrix(qc2,qh2,n2)
  m1 <- p1 * household_size$H1N1_76
  m2 <- p2 * household_size$H1N1_79
  ## adjust dimension for n=3
  na_rows <- matrix(NA, nrow = 2, ncol = ncol(m2))
  m2_adj <- rbind(m2,na_rows)
  return(cbind(m1,m2_adj))
  
}

































































