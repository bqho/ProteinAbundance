# ALL FUNCTIONS REQUIRED FOR THE ANALYSIS.
# This script should be run prior to any downstream analysis. These defined
# functions are used in subsequent scripts. 

#######-------------------------------------------------------- GetMode Function
get.mode <- function(x, y){
  x.mode = NULL
  values = log10(x[,y])

  hist.info <- hist(values,
                    breaks = seq(min(values, na.rm = T),
                                 max(values, na.rm = T),
                                 l=50+1))
  hist.mids <- hist.info$mids
  x.mode <- hist.mids[which(hist.info$counts == max(hist.info$counts))]
  
  return(10^x.mode)
}

#######---------------------------------------------- Functions for Calculations
mean.negrm <- function(x){
    values <- as.numeric(x)
    mean.value <- mean(values[which(values>0)])
    
    return(mean.value)
}

median.negrm <- function(x){
    values <- as.numeric(x)
    median.value <- median(values[which(values>0)])
    
    return(median.value)
}

var.negrm <- function(x){
    values <- as.numeric(x)
    var.value <- var(values[which(values>0)])
    
    return(var.value)
}

cv.negrm <- function(x){
    values <- as.numeric(x)
    cv.value <- cv(values[which(values>0)])
    
    return(cv.value)
}


#######-------------------------------------------------- Quantile Normalization
quant_normalize <- function(data){

# rank genes respective of their own dataset
data_rank <- apply(data, 2, rank, ties.method = "min")

# Sorting each dataset from lowest to highest, irrespective of other datasets,
# after sorting, calculate the means for each row
data_sorted <- data.frame(apply(data, 2, sort))
data_mean <- apply(data_sorted, 1, mean)

# And then we substitute the means into our ranked matrix
index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
}

data.Qnormal <- apply(data_rank, 2, index_to_mean, my_mean = data_mean)
rownames(data.Qnormal) <- rownames(data)

return(data.frame(data.Qnormal))
}

