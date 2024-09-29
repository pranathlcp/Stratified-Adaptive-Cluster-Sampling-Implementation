library(ggplot2)
library(reshape2)


###################### BEGINNING ##################### 
############ 1. Setting initial parameters ########### 

n1 = 5
n2 = 5
n = n1 + n2
N1 = 200
N2 = 200
N = N1 + N2

# Setting the condition
# Adaptive clustering will occur depending on the condition y >= condition.value
condition.value = 1

# Number of Strata
L = 2

############ 1. Setting initial parameters ########### 
######################  END ##########################



###################### BEGINNING ##################### 
############### 2. Creating the Dataset ############## 

# Creating Stratum 1 and populating the actual values
stratum.one = matrix(data = 0, nrow = 20, ncol = 10)

stratum.one[13, 4] = 2
stratum.one[14, 3] = 3
stratum.one[14, 4] = 63
stratum.one[14, 5] = 9
stratum.one[15, 4] = 16
stratum.one[15, 5] = 3
stratum.one[15, 9] = 2
stratum.one[15, 10] = 12
stratum.one[16, 9] = 2
stratum.one[16, 10] = 57
stratum.one[17, 10] = 5


# Creating Stratum 2
stratum.two = matrix(data = 0, nrow = 20, ncol = 10)

stratum.two[2, 7] = 23
stratum.two[2, 8] = 14
stratum.two[3, 7] = 36
stratum.two[3, 8] = 34
stratum.two[4, 7] = 2
stratum.two[15, 1] = 12
stratum.two[15, 2] = 1
stratum.two[16, 1] = 65
stratum.two[16, 2] = 17
stratum.two[17, 1] = 14
stratum.two[17, 2] = 5


# Creating the Population by Considering the Stratum 1 and Stratum 2
population = cbind(stratum.one, stratum.two)

## Plotting the Networks having units that satisfy the condition

# Melting for visualization purposes
pop.melted = melt(population)
pop.melted$stratum = ifelse(pop.melted$Var1 <= 10, 1, 2)
pop.melted$color = ifelse(pop.melted$value >= condition.value, "green", "#fafafa")

plot1 = ggplot(pop.melted, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = color), color = "lightgrey") +
  geom_text(data = subset(pop.melted, value >= condition.value), aes(label = value), color = "black", size = 2.5) +
  scale_fill_manual(
    name = "Legend",
    values = c("green" = "green", "#fafafa" = "#fafafa"),
    labels = c("Units not satisfying\nthe condition", "Units satisfying\nthe condition (y >= 1)")
  ) +
  theme_void() +
  geom_vline(xintercept = 10.5, color = "black", size = 1.5)+
  scale_y_reverse()

plot1

############### 2. Creating the Dataset ############## 
######################  END ##########################




#################### BEGINNING ####################### 
####### 3. Applying Stratified Random Sampling ####### 


# Random Selection (if we want, we can randomly select rather than selecting author's sample)
# stratum.one.sample = sample(length(stratum.one), n1, replace = FALSE)

# Using the author's sample
stratum.one.sample = c(27, 78, 95, 135, 149)

# Getting the sample row/col indices (coordinates)
st.one.row_indices = ((stratum.one.sample - 1) %/% ncol(stratum.one)) + 1
st.one.col_indices = ((stratum.one.sample - 1) %% ncol(stratum.one)) + 1


# Random Selection (if we want, we can randomly select rather than selecting author's sample)
# stratum.two.sample = sample(length(stratum.two), n2, replace = FALSE)

# Using the author's sample
stratum.two.sample = c(35, 44, 95, 179, 192)

# Getting the sample row/col indices (coordinates)
st.two.row_indices = ((stratum.two.sample - 1) %/% ncol(stratum.two)) + 1
st.two.col_indices = ((stratum.two.sample - 1) %% ncol(stratum.two)) + 1


# Plotting the selected samples
plot2 = plot1 + geom_rect(aes(xmin = st.one.col_indices[1] - 0.5,xmax = st.one.col_indices[1] + 0.5,ymin = st.one.row_indices[1] - 0.5,ymax = st.one.row_indices[1] + 0.5), color = "black", fill = NA, size = 1 )
plot2 = plot2 + geom_rect(aes(xmin = st.one.col_indices[2] - 0.5,xmax = st.one.col_indices[2] + 0.5,ymin = st.one.row_indices[2] - 0.5,ymax = st.one.row_indices[2] + 0.5), color = "black", fill = NA, size = 1 )
plot2 = plot2 + geom_rect(aes(xmin = st.one.col_indices[3] - 0.5,xmax = st.one.col_indices[3] + 0.5,ymin = st.one.row_indices[3] - 0.5,ymax = st.one.row_indices[3] + 0.5), color = "black", fill = NA, size = 1 )
plot2 = plot2 + geom_rect(aes(xmin = st.one.col_indices[4] - 0.5,xmax = st.one.col_indices[4] + 0.5,ymin = st.one.row_indices[4] - 0.5,ymax = st.one.row_indices[4] + 0.5), color = "black", fill = NA, size = 1 )
plot2 = plot2 + geom_rect(aes(xmin = st.one.col_indices[5] - 0.5,xmax = st.one.col_indices[5] + 0.5,ymin = st.one.row_indices[5] - 0.5,ymax = st.one.row_indices[5] + 0.5), color = "black", fill = NA, size = 1 )
plot2 = plot2 + geom_rect(aes(xmin = (st.two.col_indices+10)[1] - 0.5,xmax = (st.two.col_indices+10)[1] + 0.5,ymin = st.two.row_indices[1] - 0.5,ymax = st.two.row_indices[1] + 0.5), color = "black", fill = NA, size = 1 )
plot2 = plot2 + geom_rect(aes(xmin = (st.two.col_indices+10)[2] - 0.5,xmax = (st.two.col_indices+10)[2] + 0.5,ymin = st.two.row_indices[2] - 0.5,ymax = st.two.row_indices[2] + 0.5), color = "black", fill = NA, size = 1 )
plot2 = plot2 + geom_rect(aes(xmin = (st.two.col_indices+10)[3] - 0.5,xmax = (st.two.col_indices+10)[3] + 0.5,ymin = st.two.row_indices[3] - 0.5,ymax = st.two.row_indices[3] + 0.5), color = "black", fill = NA, size = 1 )
plot2 = plot2 + geom_rect(aes(xmin = (st.two.col_indices+10)[4] - 0.5,xmax = (st.two.col_indices+10)[4] + 0.5,ymin = st.two.row_indices[4] - 0.5,ymax = st.two.row_indices[4] + 0.5), color = "black", fill = NA, size = 1 )
plot2 = plot2 + geom_rect(aes(xmin = (st.two.col_indices+10)[5] - 0.5,xmax = (st.two.col_indices+10)[5] + 0.5,ymin = st.two.row_indices[5] - 0.5,ymax = st.two.row_indices[5] + 0.5), color = "black", fill = NA, size = 1 )

plot2

####### 3. Applying Stratified Random Sampling ####### 
#######################  END #########################




################### BEGINNING ####################### 
####### 4. Applying Adaptive Cluster Sampling ####### 

st.one.sampled_indices = cbind(row = st.one.row_indices, col = st.one.col_indices)
st.two.sampled_indices = cbind(row = st.two.row_indices, col = st.two.col_indices + 10)
pop.sampled_indices_df = as.data.frame(rbind(st.one.sampled_indices, st.two.sampled_indices))

# Function to get neighborhood indices
get_neighborhood = function(row, col) {
  neighborhood = data.frame(
    row = c(row - 1, row + 1, row, row),
    col = c(col, col, col - 1, col + 1)
  )
  return(neighborhood)
}

added_new_cells = TRUE
while (added_new_cells) {
  added_new_cells = FALSE
  
  # Iterating through each cell in pop.sampled_indices_df
  for (i in 1:nrow(pop.sampled_indices_df)) {
    cell_row = pop.sampled_indices_df$row[i]
    cell_col = pop.sampled_indices_df$col[i]
    
    # Getting the neighborhood for the current cell
    neighborhood = get_neighborhood(cell_row, cell_col)
    
    # Filtering neighborhood cells that fall within the bounds of the population matrix
    neighborhood = neighborhood[neighborhood$row >= 1 & neighborhood$row <= nrow(population) &
                                  neighborhood$col >= 1 & neighborhood$col <= ncol(population), ]
    
    # Checking each cell in the filtered neighborhood
    for (j in 1:nrow(neighborhood)) {
      neighbor_row = neighborhood$row[j]
      neighbor_col = neighborhood$col[j]
      
      if (!(population[neighbor_row, neighbor_col] >= condition.value)) {
        
        if (population[cell_row,cell_col] >= condition.value) {
          pop.sampled_indices_df = rbind(pop.sampled_indices_df, data.frame(row = neighbor_row, col = neighbor_col))
        }
        
        next
      }
      
      if (population[neighbor_row, neighbor_col] >= condition.value) {
        pop.sampled_indices_df = rbind(pop.sampled_indices_df, data.frame(row = neighbor_row, col = neighbor_col))
      }
      
      # Getting the neighborhood for the current cell
      neighbors.neighborhood = get_neighborhood(neighbor_row, neighbor_col)
      
      neighbors.neighborhood = neighbors.neighborhood[
        neighbors.neighborhood$row >= 1 & 
          neighbors.neighborhood$row <= nrow(population) &
          neighbors.neighborhood$col >= 1 & 
          neighbors.neighborhood$col <= ncol(population), 
      ]
      
      for (k in 1:nrow(neighbors.neighborhood)) {
        neighbors.neighbor_row = neighbors.neighborhood$row[k]
        neighbors.neighbor_col = neighbors.neighborhood$col[k]
        
        # Checking if the neighbor is not already in sampled indices and satisfies the condition
        if (!any(pop.sampled_indices_df$row == neighbors.neighbor_row &
                 pop.sampled_indices_df$col == neighbors.neighbor_col) &&
            population[neighbors.neighbor_row, neighbors.neighbor_col] >= condition.value) {
          
          # Adding the new cell to pop.sampled_indices_df
          pop.sampled_indices_df = rbind(pop.sampled_indices_df, data.frame(row = neighbors.neighbor_row, col = neighbors.neighbor_col))
          
          # Updating the flag to continue the loop
          added_new_cells = TRUE
        }
      }
    }
  }
}

# Removing the duplicates
pop.sampled_indices_df = unique(pop.sampled_indices_df)
pop.sampled_indices_df
dim(pop.sampled_indices_df)


# Plotting the final sample
pop.melted$highlight = ifelse(
  interaction(pop.melted$Var1, pop.melted$Var2) %in% interaction(pop.sampled_indices_df$row, pop.sampled_indices_df$col),
  "highlight",
  ifelse(pop.melted$value >= condition.value, "network", "other")
)

plot3 = ggplot(pop.melted, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = highlight), color = "lightgrey") +
  geom_text(data = subset(pop.melted, value >= condition.value), aes(label = value), color = "black", size = 2.5) +
  scale_fill_manual(
    name = "Legend",
    values = c("green" = "green", "#fafafa" = "#fafafa", "highlight" = "red", "network" = "green", "other" = "#fafafa"),  # Include red for highlighted cells
    labels = c("Units added to the\nfinal sample", "Units satisfying the condition\nbut not in the final sample", "Units not selected\nto the final sample")
  ) +
  theme_void() +
  geom_vline(xintercept = 10.5, color = "black", size = 1.5) +
  scale_y_reverse()

plot3


# Adding the values to the pop.sampled_indices_df
pop.sampled_indices_df$value = population[cbind(pop.sampled_indices_df$row, pop.sampled_indices_df$col)]

# Ordering the pop.sampled_indices_df
pop.sampled_indices_df = pop.sampled_indices_df[with(pop.sampled_indices_df, order(row, col)), ]

####### 4. Applying Adaptive Cluster Sampling #######
#######################  END #########################





################### BEGINNING ####################### 
################# 5. Data Processing ################ 

# Populating the Network Column
# Ideally, this should be coded such that network numbers can be obtained dynamically
pop.sampled_indices_df$network = 0

network_1_coords = cbind(c(13, 14, 14, 14, 15, 15), c(4, 4, 5, 3, 4, 5))
network_2_coords = cbind(
  c(15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17),
  c(9, 10, 11, 12, 9, 10, 11, 12, 10, 11, 12)
)

pop.sampled_indices_df$network[rowSums(outer(pop.sampled_indices_df$row, network_1_coords[,1], "==") &
                                         outer(pop.sampled_indices_df$col, network_1_coords[,2], "==")) > 0] = 1

pop.sampled_indices_df$network[rowSums(outer(pop.sampled_indices_df$row, network_2_coords[,1], "==") &
                                         outer(pop.sampled_indices_df$col, network_2_coords[,2], "==")) > 0] = 2

# Adding the Stratum Column to pop.sampled_indices_df
pop.sampled_indices_df$stratum = ifelse(pop.sampled_indices_df$col <= 10, 1, 2)

# Adding an indicator to determine whether the unit was originally/initially sampled
pop.sampled_indices_df$initially_sampled = FALSE
for(i in 1:nrow(pop.sampled_indices_df)) {
  final_sample_row = pop.sampled_indices_df[i, "row"]
  final_sample_col = pop.sampled_indices_df[i, "col"]
  
  for (j in 1:nrow(st.one.sampled_indices)) {
    initial_sample_row = st.one.sampled_indices[j, "row"]
    initial_sample_col = st.one.sampled_indices[j, "col"]
    if (final_sample_row == initial_sample_row && final_sample_col == initial_sample_col) {
      pop.sampled_indices_df[i, "initially_sampled"] = TRUE
    }
  }
  
  for (j in 1:nrow(st.two.sampled_indices)) {
    initial_sample_row = st.two.sampled_indices[j, "row"]
    initial_sample_col = st.two.sampled_indices[j, "col"]
    if (final_sample_row == initial_sample_row && final_sample_col == initial_sample_col) {
      pop.sampled_indices_df[i, "initially_sampled"] = TRUE
    }
  }
}

# Processed samples df
pop.sampled_indices_df

################# 5. Data Processing ################ 
#######################  END ########################



################## BEGINNING #######################
######## 6. Computing the estimator t1 #############
####### Conventional Stratified Random Sampling ####

n1.sample = pop.sampled_indices_df[pop.sampled_indices_df$initially_sampled &
                                     pop.sampled_indices_df$stratum == 1, ]
n2.sample = pop.sampled_indices_df[pop.sampled_indices_df$initially_sampled &
                                     pop.sampled_indices_df$stratum == 2, ]



# sample size, sample mean, sample variance of each strata
strata.N = Wh = strata.n = strata.avg = strata.var = rep(NA, L)

for(i in 1:L) {
  strata.N[i] = N1
  strata.n[i] = n1
  Wh[i] = N1/N
  strata.sample.SRSWOR = pop.sampled_indices_df[pop.sampled_indices_df$initially_sampled &
                                                  pop.sampled_indices_df$stratum == i, ]

  strata.avg[i] = mean(strata.sample.SRSWOR$value, na.rm = TRUE)
  strata.var[i] = var(strata.sample.SRSWOR$value, na.rm = TRUE)
}


result.stratified = data.frame(Stratum = c(1,2), Nh = strata.N, nh = strata.n, Average = strata.avg, Variance = strata.var, Wh = Wh)
result.stratified

# t1 = unbiased mean estimator = 1.1
t1 = sum(result.stratified$Wh*result.stratified$Average)
t1

# var.t1 = unbiased variance estimator for y_bar = 0.741
var.t1 = sum( (result.stratified$Wh)^2 * ( 1 - result.stratified$nh/result.stratified$Nh) * result.stratified$Variance/result.stratified$nh )
var.t1

########### 6. Computing the estimator t1 ###############
######################### END ###########################





##################### BEGINNING ########################
########### 7. Computing the estimator t2 ##############

strata.N = strata.n = strata.avg = strata.w_h = strata.var = rep(NA, L)

for(h in 1:L) {
  stratum.data = pop.sampled_indices_df[which(pop.sampled_indices_df$stratum == h), ]
  
  w_h = 0
  sampled_nh = pop.sampled_indices_df[pop.sampled_indices_df$initially_sampled == TRUE 
                                      & pop.sampled_indices_df$stratum == h, ]
  
  strata.N[h] = dim(pop.melted[pop.melted$stratum == h, ])[1]
  strata.n[h] = dim(sampled_nh)[1]
  
  wh_list = rep(0, strata.n[h])
  
  for(i in 1:nrow(sampled_nh)) {
    sampled_unit = sampled_nh[i, ]
    if (sampled_unit$network == 0) {
      
      # By definition
      w_h = w_h + 0
      wh_list[h] = 0
      next
    }
    
    networked_data = pop.sampled_indices_df[which(pop.sampled_indices_df$network == sampled_unit$network), ]
    w_hi = sum(networked_data$value)/nrow(networked_data)
    w_h = w_h + w_hi
    wh_list[i] = w_hi
  }
  
  strata.w_h[h] = w_h
  strata.var[h] = var(wh_list)
}

result = data.frame(Nh = strata.N, nh = strata.n, w_h =  strata.w_h, strata.var = strata.var)
result

# t2 = mu1.hat = 3.345455
t2 = (1/N) * sum( (result$Nh/result$nh) * result$w_h)
t2

# var_t2 = var_mu1.hat = 4.104992
var_t2 = (1/N^2) * sum(result$Nh*(result$Nh - result$nh) * (result$strata.var/result$nh) )
var_t2

########### 7. Computing the estimator t2 ##############
######################## END ###########################




#################### BEGINNING ##########################
###### 8. Computing the estimator t2.doublePrime ########

strata.N = strata.n = strata.avg = strata.w_h_doublePrime = strata.var = rep(NA, L)

for(i in 1:L) {
  stratum.data = pop.sampled_indices_df[which(pop.sampled_indices_df$stratum == i), ]
  
  w_h_doublePrime = 0
  sampled_nh = pop.sampled_indices_df[pop.sampled_indices_df$initially_sampled == TRUE 
                                      & pop.sampled_indices_df$stratum == i, ]
  
  strata.N[i] = dim(pop.melted[pop.melted$stratum == i, ])[1]
  strata.n[i] = dim(sampled_nh)[1]
  
  wh_DoublePrimes = rep(0, strata.n[i])
  
  for(j in 1:nrow(sampled_nh)) {
    sampled_unit = sampled_nh[j, ]
    if (sampled_unit$network == 0) {
      
      # By definition
      w_h_doublePrime = w_h_doublePrime + 0
      wh_DoublePrimes[i] = 0
      next
    }
    
    intersected_data = stratum.data[which(stratum.data$network == sampled_unit$network), ]
    computation = sum(intersected_data$value)/nrow(intersected_data)
    w_h_doublePrime = w_h_doublePrime + computation
    wh_DoublePrimes[j] = computation
  }
  
  strata.w_h_doublePrime[i] = w_h_doublePrime
  strata.var[i] = var(wh_DoublePrimes)
}

result = data.frame(Nh = strata.N, nh = strata.n, w_h_doublePrime =  strata.w_h_doublePrime, strata.var = strata.var)
result

# t2.doublePrime = mu1.hat.doublePrime = 3.16
t2.doublePrime = sum(result$w_h_doublePrime*result$Nh/result$nh)/N
t2.doublePrime

# var_t2.doublePrime = 3.65196
var_t2.doublePrime = (1/N^2) * sum(result$Nh*(result$Nh - result$nh) * (result$strata.var/result$nh) )
var_t2.doublePrime


###### 8. Computing the estimator t2.doublePrime ########
######################### END ###########################







##################### BEGINNING ##########################
############# 9. Computing the estimator t3 ##############

# We have two networks having units that satisfy our condition (y >= 1)
# As a result, we can compute two inclusion probabilities.
# For units that don't satisfy the condition, the inclusion probability will be pi.0

# pi.0 = 0.025 
# Same for each stratum because of equal sample and stratum sizes
pi.0 = n1/N1
pi.0

# A function for computing the inclusion probability
get_inclusion_probability = function(network) {
  product_computation = 1
  for (k in 1:L) {
    N.k = dim(pop.melted[pop.melted$stratum == k, ])[1]
    n.k = dim(pop.sampled_indices_df[pop.sampled_indices_df$stratum == k &
                                       pop.sampled_indices_df$initially_sampled, ])[1]
    
    network_stratum_intersection = pop.sampled_indices_df[pop.sampled_indices_df$stratum == k &
                                                            pop.sampled_indices_df$network == network, ]
    
    x.ki = dim(network_stratum_intersection)[1]
    
    computation = choose(N.k - x.ki, n.k)/choose(N.k, n.k)
    
    product_computation = product_computation * computation
  }
  pi = 1 - product_computation
}

# A function for computing the inclusion probability
get_joint_inclusion_probability = function(networks, pi.i, pi.j) {
  i = networks[1]
  j = networks[2]
  product_computation = 1
  for (k in 1:L) {
    N.k = dim(pop.melted[pop.melted$stratum == k, ])[1]
    n.k = dim(pop.sampled_indices_df[pop.sampled_indices_df$stratum == k &
                                       pop.sampled_indices_df$initially_sampled, ])[1]
    
    ki_intersection = pop.sampled_indices_df[pop.sampled_indices_df$stratum == k &
                                                            pop.sampled_indices_df$network == i, ]
    kj_intersection = pop.sampled_indices_df[pop.sampled_indices_df$stratum == k &
                                               pop.sampled_indices_df$network == j, ]
    
    x.ki = dim(ki_intersection)[1]
    x.kj = dim(kj_intersection)[1]
    
    computation = choose(N.k - x.ki - x.kj, n.k)/choose(N.k, n.k)
    
    product_computation = product_computation * computation
  }
  
  p.ij = 1 - (1 - pi.i) - (1 - pi.j) + product_computation
}

# pi.1 = 0.1426134
# pi.2 = 0.2455432
# pi.12 = 0.03240313
pi.1 = get_inclusion_probability(network = 1)
pi.2 = get_inclusion_probability(network = 2)
pi.1
pi.2

pi.12 = get_joint_inclusion_probability(c(1,2), pi.1, pi.2)
pi.12

## t3 Computation

# Creating z.i indicator variable
pop.sampled_indices_df$z = ifelse(pop.sampled_indices_df$initially_sampled & pop.sampled_indices_df$network != 0, 1, 0)
pop.sampled_indices_df$pi = ifelse(pop.sampled_indices_df$network == 1, pi.1, 
                                   ifelse(pop.sampled_indices_df$network == 2, pi.2, pi.0))

summation = 0
for (i in 1:nrow(pop.sampled_indices_df)) {
  unit = pop.sampled_indices_df[i, ]
  y_value = sum(pop.sampled_indices_df[pop.sampled_indices_df$network == unit$network, "value"])
  summation = summation + (y_value * unit$z) / unit$pi
}

# t3 = 3.63772
t3 = (1/N) * summation
t3

## var(t3) Computation
outer_sum = 0
for (i in 1:nrow(pop.sampled_indices_df)) {
  inner_sum = 0
  unit.i = pop.sampled_indices_df[i, ]
  y.i = sum(pop.sampled_indices_df[pop.sampled_indices_df$network == unit.i$network, "value"])
  
  for (j in 1:nrow(pop.sampled_indices_df)) {
    unit.j = pop.sampled_indices_df[j, ]
    y.j = sum(pop.sampled_indices_df[pop.sampled_indices_df$network == unit.j$network, "value"])
    pi.ij = pi.12
    
    if (i == j) {
      pi.ij = unit.i$pi
    }
    
    comp = ((y.i * y.j * unit.i$z * unit.j$z) / pi.ij) * ( (pi.ij / (unit.i$pi * unit.j$pi) ) - 1 )
    
    inner_sum = inner_sum + comp
  }
  outer_sum = outer_sum + inner_sum
}

# var.t3 = 4.780367
var.t3 = (1/N^2) * outer_sum
var.t3


############# 9. Computing the estimator t3 ##############
########################### END ##########################



