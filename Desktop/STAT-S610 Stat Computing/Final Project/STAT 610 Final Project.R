
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
  D_star1[is.na(D_star1)] <- 0
  D_star2[is.na(D_star2)] <- 0
  D1[is.na(D1)] <- 0
  D2[is.na(D2)] <- 0
  norm1 <- norm(D1 - D_star1, type = "F")
  norm2 <- norm(D2 - D_star2, type = "F")
  dist <- 0.5 * (norm1 + norm2)
  return(dist) 
  
}

household_size <- list(A_78 <- colSums(H3N2_78,na.rm = T),
                       A_81 <-colSums(H3N2_81,na.rm = T),
                       B_76 <- colSums(H1N1_76,na.rm = T),
                       B_79 <- colSums(H1N1_79,na.rm = T))
names(household_size) <- c('A_78','A_81','B_76','B_79')

calculate_w_matrix <- function(qc, qh, n) {
  # Initialize the matrix with NA values
  w_matrix <- matrix(0, nrow = n+1, ncol = n)
  
  # Set the first row (w_0s)
  for (s in 1:n){
    w_matrix[1,s] <- qc^(s)
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
  m1 <- sweep(p1, 2, household_size$A_78, `*`)
  m2 <- sweep(p2, 2, household_size$A_81, `*`)
  return(cbind(m1,m2))
}

#Initialize the prior and 
N = 5e+5
qc1 = runif(N)
qh1 = runif(N)
qc2 = runif(N)
qh2 = runif(N)
D1 = H3N2_78 
D2 = H3N2_81

#simulated_data <- make_table_2(0.6,0.7,0.4,0.5)
#Table 2
generate_abc <- function(qc1, qc2, qh1, qh2, N, epsilon, D1, D2){
  accepted_sample = data.frame(qc1 = numeric(0),
                               qh1 = numeric(0),
                               qc2 = numeric(0),
                               qh2 = numeric(0))
  
  for (i in 1:N){
    #simulate data for given parameters
    simulated_data <- make_table_2(qc1[i], qh1[i], qc2[i], qh2[i])
    #assign the simulated data to a new variable 
    D_star1 <- simulated_data[,1:5]
    D_star2 <- simulated_data[,6:10]
    #calculate the distance
    dist <- distance(D1, D2, D_star1, D_star2)
    
    #Accept/Reject Condition
    if (dist < epsilon){
      new_values <- c(qc1[i], qh1[i], qc2[i], qh2[i])
      accepted_sample <- rbind(accepted_sample, new_values)
    }
  }
  colnames(accepted_sample) <- c("qc1", "qh1", "qc2", "qh2")
  return(accepted_sample)
  
}

accepted_sample <- generate_abc(qc1,qc2,qh1,qh2,N,15, D1, D2)
#accepted_sample

long_df <- pivot_longer(accepted_sample, 
                        cols = c(qc1, qh1, qc2, qh2), 
                        names_to = c(".value", "group"), 
                        names_pattern = "(..)(.)")
#long_df

ggplot(long_df, aes(x = qh, y = qc, color = group)) +
  geom_point() +
  xlim(0,1)+ylim(0,1)+
  labs(title = "Scatter Plot by Group",
       x = "Qh Axis",
       y = "Qc Axis") +
  theme_minimal()


# Table 3
make_table_3 <- function(qc1,qh1,qc2,qh2){
  ## for table 3, 75-76 has col 5 and 78-79 has col 3
  n1 <- 5
  n2 <- 3
  p1 <- calculate_w_matrix(qc1,qh1,n1)
  p2 <- calculate_w_matrix(qc2,qh2,n2)
  m1 <- sweep(p1, 2, household_size$B_76, `*`)
  m2 <- sweep(p2, 2, household_size$B_79, `*`)
  ## adjust dimension for n=3
  na_rows <- matrix(NA, nrow = 2, ncol = ncol(m2))
  m2_adj <- rbind(m2,na_rows)
  return(cbind(m1,m2_adj))

}
#Initialize the prior and
N = 3e+5
qc1 = runif(N)
qh1 = runif(N)
qc2 = runif(N)
qh2 = runif(N)
D3 = H1N1_76
D4 = H1N1_79

generate_abc2 <- function(qc1, qc2, qh1, qh2, N, epsilon, D1, D2){
  accepted_sample = data.frame(qc1 = numeric(0),
                               qh1 = numeric(0),
                               qc2 = numeric(0),
                               qh2 = numeric(0))

  for (i in 1:N){
    #simulate data for given parameters
    simulated_data <- make_table_3(qc1[i], qh1[i], qc2[i], qh2[i])
    #assign the simulated data to a new variable
    D_star1 <- simulated_data[,1:5]
    D_star2 <- simulated_data[1:4,6:8]
    #calculate the distance
    dist <- distance(D1, D2, D_star1, D_star2)

    #Accept/Reject Condition
    if (dist < epsilon){
      new_values <- c(qc1[i], qh1[i], qc2[i], qh2[i])
      accepted_sample <- rbind(accepted_sample, new_values)
    }
  }
  colnames(accepted_sample) <- c("qc1", "qh1", "qc2", "qh2")
  return(accepted_sample)

}

accepted_sample2 <- generate_abc2(qc1,qc2,qh1,qh2,N,7, D3, D4)
accepted_sample2

long_df2 <- pivot_longer(accepted_sample2,
                         cols = c(qc1, qh1, qc2, qh2),
                         names_to = c(".value", "group"),
                         names_pattern = "(..)(.)")
long_df2

ggplot(long_df2, aes(x = qh, y = qc, color = group)) +
  geom_point() + xlim(0,1) + ylim(0,1) +
  labs(title = "Scatter Plot by Group",
       x = "Qh Axis",
       y = "Qc Axis") +
  theme_minimal() +
  geom_encircle(aes(group = group), 
                expand = unit(0.1, "cm"))











library(testthat)
# Test for distance function
test_that("distance function outputs correct Frobenius norm", {
  # Create two known matrices
  D1 <- matrix(c(1, 2, 3, 4), nrow = 2)
  D2 <- matrix(c(5, 6, 7, 8), nrow = 2)
  D_star1 <- matrix(c(1, 1, 1, 1), nrow = 2)
  D_star2 <- matrix(c(2, 2, 2, 2), nrow = 2)
  
  # Manually calculate Frobenius norm
  manual_norm <- (sqrt(sum((D1 - D_star1)^2)) + sqrt(sum((D2 - D_star2)^2))) / 2
  
  # Check if the function's output matches manual calculation
  expect_equal(distance(D1, D2, D_star1, D_star2), manual_norm)
})

# Test for calculate_w_matrix function
test_that("calculate_w_matrix outputs a matrix with correct properties", {
  n <- 5
  qc <- 0.3
  qh <- 0.7
  
  # Obtain the matrix from the function
  w_matrix <- calculate_w_matrix(qc, qh, n)
  
  # Check that all column sums to one
  expect_equal(unname(colSums(w_matrix)), rep(1, n))
  
  # Check that each entry is in [0, 1]
  expect_true(all(w_matrix >= 0 & w_matrix <= 1))
})

# Test for make_table_2 function
test_that("make_table_2 outputs a table with weighted column sums", {
  qc1 <- 0.2
  qh1 <- 0.8
  qc2 <- 0.4
  qh2 <- 0.6
  # Manually create weights and calculate expected table
  expected_table <- matrix(c(15.8, 4.20,  0.384,  0.0704, 0.00352, 21.6, 13.440,  
                             3.84000, 1.587200,  0.1945600, 63.2, 26.88,  2.94912,
                             0.5767168, 0.02883584, 32.4, 24.192,  6.2208,  2.057011,
                             0.1891123, 0.0, 73.92, 12.97613,  3.0450647, 0.16240345,
                             0.0, 46.368, 14.30784,  4.258013,  0.31317, 0.0, 0.00, 
                             31.69075, 11.8988210, 0.76152455,  0.0,  0.00, 35.63136,
                             12.724671,  0.8422911, 0.0,  0.0,  0.0, 28.4089975, 2.9090,
                             0.0,  0.0,  0.0, 41.373104,  3.2863591, 0.0,  0.0,  0.0,  
                             0.0, 7.13463482, 0.0,  0.0,  0.0,  0.0, 14.1745), 
                           nrow = 6,ncol = 10, byrow = T)  
  # Obtain the table from the function
  result_table <- unname(make_table_2(qc1, qh1, qc2, qh2))
  
  # Check if the function's output matches the expected table
  tolerance <- 1e-5  # You can adjust the tolerance level as needed
  expect_equal(result_table, expected_table, tolerance = tolerance)
})































