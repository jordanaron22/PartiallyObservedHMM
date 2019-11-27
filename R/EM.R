##################################################
############## Simulation Functions ##############
##################################################
#' Generates Simulated Latent Data
#'
#' Generates simulated latent data from initial and transition probabilities.  Incorporated in GenerateSimulatedData()
#'
#' @param n number of individuals to generate data for
#' @param t number of visits to generate data for
#' @param p probability of being a stayer in the healthy state for latent data
#'
#' @return list of true data and vector of indicators for whether each individual is a stayer
#'
#' @export
GenerateData <- function(n,t,p){
  init_hc <- c(0.272,0.000,0.116,0.070,0.003,0.001,0.047,0.025,0.106,0.057,0.014,0.020,0.047,0.025,0.106,0.057,0.014,0.020)
  tran_hc <- matrix(c(0.809,0.000,0.089,0.000,0.002,0.000,0.073,0.000,0.026,0.000,0.001,0.0,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.000,0.130,0.000,0.309,0.000,0.000,0.000,0.426,0.000,0.136,0.000,0.0,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.643,0.000,0.260,0.000,0.007,0.000,0.043,0.000,0.045,0.000,0.002,0.0,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.000,0.000,0.000,0.751,0.000,0.002,0.000,0.128,0.000,0.119,0.000,0.0,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.554,0.000,0.290,0.000,0.071,0.000,0.031,0.000,0.054,0.000,0.000,0.0,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.000,0.000,0.000,0.500,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.5,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.596,0.000,0.056,0.000,0.001,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.206,0.000,0.131,0.000,0.010,0.000,
                      0.000,0.000,0.000,0.164,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.509,0.000,0.304,0.000,0.023,
                      0.367,0.000,0.103,0.000,0.011,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.167,0.000,0.321,0.000,0.031,0.000,
                      0.000,0.000,0.000,0.160,0.000,0.002,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.274,0.000,0.499,0.000,0.065,
                      0.357,0.000,0.152,0.000,0.015,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.134,0.000,0.231,0.000,0.111,0.000,
                      0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.316,0.000,0.474,0.000,0.211,
                      0.435,0.000,0.043,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.320,0.003,0.186,0.002,0.011,0.000,
                      0.000,0.000,0.000,0.089,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.560,0.000,0.325,0.000,0.025,
                      0.285,0.000,0.080,0.000,0.005,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.187,0.002,0.395,0.003,0.042,0.000,
                      0.000,0.000,0.000,0.100,0.000,0.003,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.246,0.000,0.589,0.000,0.063,
                      0.207,0.000,0.136,0.000,0.007,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.121,0.001,0.405,0.003,0.118,0.001,
                      0.000,0.000,0.000,0.042,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.042,0.000,0.583,0.000,0.333),18,18, byrow = T)


  data_mat <- matrix(0L, nrow = n, ncol = t)
  stayer_vec <- stats::rbinom(n,1,p)

  for (i in 1:n){
    if (stayer_vec[i] == 0){
      data_mat[i,1] <- which(rmultinom(1,1,init_hc) == 1) - 1
      for (j in 2:t){
        p_current <- tran_hc[data_mat[i,j-1]+1,]
        data_mat[i,j] <- which(rmultinom(1,1,p_current) == 1) - 1
      }
    }
  }


  data_true_array <- array(NaN, dim = c(dim(data_mat)[1], dim(data_mat)[2], 2))
  data_true_array[,,1] <- data_mat
  return(list(data_mat, stayer_vec))
}

#' Generates Simulated Observed Data
#'
#' Generates simulated observed data from simulated latent data by adding measurment error.  Incorporated in GenerateSimulatedData()
#'
#' @param data_true simulated latent data that error is added to
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#'
#' @return list of observed data and true classification matrix
#'
#'
#' @export
GenerateObs <- function(data_true,max2,max3){
  class_hc <- GetClass()
  data_obs <- matrix(NA,dim(data_true)[1], dim(data_true)[2])
  k1 <- length(numOfStates(max2,max3))
  class_true <- matrix(0, k1,k1)
  for (i in 1:dim(data_true)[1]){
    to_treatment <- F
    for (j in 1:dim(data_true)[2]){
      p_current <- class_hc[data_true[i,j]+1,]
      new_val <- which(rmultinom(1,1,p_current) == 1) - 1
      if (!to_treatment){
        class_true[data_true[i,j]+1, new_val+1] <- class_true[data_true[i,j]+1, new_val+1] + 1
        data_obs[i,j] <- new_val
      }

      if (new_val %% 2 == 1){
        to_treatment <- T
      }

    }
  }
  return(list(data_obs, class_true))
}

#' Seperates Three Tests
#'
#' Takes in data of simulated latent data of three tests and seperates them into three seperate matrices.  Incorporated in GenerateSimulatedData()
#'
#' @param data_obs simulated observed latent data
#'
#' @return list of three data matrices
#'
#' @export
Seperator <- function(data_obs){
  hpv <- matrix(NA, dim(data_obs)[1], dim(data_obs)[2])
  cyt <- matrix(NA, dim(data_obs)[1], dim(data_obs)[2])
  colpo <- matrix(NA, dim(data_obs)[1], dim(data_obs)[2])

  for (i in 1:dim(data_obs)[1]){
    for (j in 1:dim(data_obs)[2]){
      if (!is.na(data_obs[i,j])){
        hpv[i,j] <- data_obs[i,j] %/% 6
        cyt[i,j] <- (data_obs[i,j] %% 6) %/% 2
        colpo[i,j] <- data_obs[i,j] %% 2
      }
    }
  }
  return(list(hpv,cyt,colpo))
}

#' Introduces Misisng data
#'
#' Introduces missingness that is independent of the value of the other two tests.  Incorporated in GenerateSimulatedData()
#'
#' @param ind_data_obs data to introduce independent randomness to
#' @param p probability that each observation is missing
#'
#' @return data matrix with missing values
#'
#' @export
IntroduceMissing <- function(ind_data_obs, p){
  for (i in 1:dim(ind_data_obs)[1]){
    for (j in 1:dim(ind_data_obs)[2]){
      if (rbinom(1,1,p)){
        ind_data_obs[i,j] <- NA
      }
    }
  }
  return(ind_data_obs)
}

#' Introduces Misisng data
#'
#' Introduces missingness that is dependent on the value of the other two tests.  Incorporated in GenerateSimulatedData()
#'
#' @param mat1 first dependent test
#' @param mat2 second dependent test
#' @param mat3 test that misssingness is added to
#'
#' @return data matrix with missing values
#'
#' @export
IntroduceMissingColpo <- function(mat1, mat2, mat3){
  p_matrix <- matrix(c(0.845,0.833,0.071,0.500,
                       0.889,0.840,0.091,1.000,
                       0.843,0.778,0.128,1.000,
                       0.664,0.518,0.091,0.375,
                       0.820,0.745,0.083,0.999), 5,4,byrow = T)
  for (i in 1:dim(mat3)[1]){
    for (j in 1:dim(mat3)[2]){

      if (!is.na(mat1[i,j])){
        mat1_ind <- mat1[i,j] + 1
      } else {
        mat1_ind <- 5
      }

      if (!is.na(mat2[i,j])){
        mat2_ind <- mat2[i,j] + 1
      } else {
        mat2_ind <- 4
      }

      if (rbinom(1,1,p_matrix[mat1_ind,mat2_ind])){
        mat3[i,j] <- NA
      }

    }
  }
  return(mat3)
}

#' Introduced Missing Data
#'
#' Introduces missing data to all three tests at the same time point to mimic mixed visits. Incorporated in GenerateSimulatedData()
#'
#' @param mat1 test 1
#' @param mat2 test 2
#' @param mat3 test 2
#' @param p_time vector of length t for missing probabilities for each time.
#'
#' @return list of three matrices with concurrent missing data
#'
#' @export
IntroduceMissingAll <- function(mat1,mat2,mat3, p_time){
  for (i in 1:dim(mat1)[1]){
    for (j in 1:dim(mat1)[2]){
      if (rbinom(1,1,p_time[j])){
        mat1[i,j,] <- NA
        mat2[i,j] <- NA
        mat3[i,j] <- NA
      }
    }
  }
  return(list(mat1,mat2,mat3))
}

#' Calculates True Intial Probabilities
#'
#' Calculates True Intial Probabilities from simulated latent data.  Incorporated in GenerateSimulatedData()
#'
#' @param data_true simulated latent data
#' @param stayer_vec_true vector of latent indicators for whether an individual is a stayer
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#'
#' @return true initial vector
#'
#' @export
GetInitTrue <- function(data_true,stayer_vec_true,max2,max3){
  init_counts_true <- numeric(length(numOfStates(max2,max3)))
  for (i in 1:dim(data_true)[1]){
    if (stayer_vec_true[i] == 0 & !is.na(data_true[i,1])){
      init_counts_true[data_true[i,1] + 1] <- init_counts_true[data_true[i,1] + 1] + 1
    }
  }
  init_true <- init_counts_true/sum(init_counts_true)
  return(init_true)
}

#' Calculates True Transiition Probabilities
#'
#' Calculates True Transiition Probabilities from simulated latent data.  Incorporated in GenerateSimulatedData()
#'
#' @param data_true simulated latent data
#' @param stayer_vec_true vector of latent indicators for whether an individual is a stayer
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#'
#' @return true transition matrix
#'
#' @export
GetTranTrue <-  function(data_true,stayer_vec_true,max2,max3){
  tran_counts_true <- matrix(0,length(numOfStates(max2,max3)),length(numOfStates(max2,max3)))
  for (i in 1:dim(data_true)[1]){
    for (j in 2:dim(data_true)[2]){
      if (stayer_vec_true[i] == 0 & !is.na(data_true[i,j-1]) & !is.na(data_true[i,j])){
        tran_counts_true[data_true[i,j-1] + 1, data_true[i,j] + 1] <- tran_counts_true[data_true[i,j-1] + 1, data_true[i,j] + 1] + 1
      }
    }
  }

  tran_true <- Normalize(tran_counts_true)
  return(tran_true)
}

#' Introduces Persistence
#'
#' Introduces one observation persistence to binary test.  Adds time dependence while still maintaining Markov identity.  Incorporated in GenerateSimulatedData.
#'
#' @param hpv binary data matrix where one-obersvation back persistence is introduced
#'
#' @return time dependent trinary data matrix
#'
#' @export
IntroducePersistence <- function(hpv){
  for (i in 1:dim(hpv)[1]){
    for (j in 1:dim(hpv)[2]){
      if (!is.na(hpv[i,j])){
        if (j == 1){hpv[i,j] <- persitenceHelper(hpv[i,j], NA)}
        else {hpv[i,j] <- persitenceHelper(hpv[i,j], hpv[i,j-1])}
      }
    }
  }

  hpv_pers <- array(NA, dim = c(dim(hpv)[1], dim(hpv)[2], 3))
  for (i in 1:dim(hpv)[1]){
    for (j in 1:dim(hpv)[2]){
      if (!is.na(hpv[i,j])){
        if (hpv[i,j] < 3){hpv_pers[i,j,1] <- hpv[i,j]}
        else if (hpv[i,j] == 3){
          hpv_pers[i,j,1] <- 1
          hpv_pers[i,j,2] <- 2
        }
      }
    }
  }
  return(hpv_pers)
}

#' Helper Function for IntroducePersistence
#'
#' Helper Function for IntroducePersistence
#'
#' @param val current value that is being made time dependent
#' @param one_back value from last time point
#'
#' @return time dependent value
#'
#' @export
persitenceHelper <- function(val, one_back){
  if (val == 0){return(0)}

  else if (!is.na(one_back)){
    if (one_back == 0){return(1)}
    else if (one_back == 1){return(2)}
    else if (one_back == 2){return(2)}
    else if (one_back == 3){return(2)}
  }
  else{return(3)}

}

#' Generates Simulated Data
#'
#' Generates simulated HPV, cytology and colposcopy data by using the following functions:
#'    GenerateData, GenerateObs, GetClass, GetInitTrue, GetTranTrue, IntroduceMissing, IntroduceMissingAll, IntroduceMissingColpo, and Seperator
#'
#' @param n number of individuals to generate
#' @param t number of timepoints to generate
#' @param p probability of being a stayer for simulated data
#' @param max2 number of possibilities from test 2
#' @param max3 number of possibilities from test
#'
#' @return list of three matrices, one for each test
#'
#' @export
GenerateSimulatedData <- function(n,t,p,max2,max3){
  data_true_full <- GenerateData(n,t,p)
  data_true <- data_true_full[[1]]
  stayer_vec_true <- data_true_full[[2]]

  data_obs_full <- GenerateObs(data_true,max2,max3)
  data_obs <- data_obs_full[[1]]
  class_counts_true <- data_obs_full[[2]]


  init_true <- GetInitTrue(data_true,stayer_vec_true,max2,max3)
  tran_true <- GetTranTrue(data_true,stayer_vec_true,max2,max3)
  class_true <- Normalize(class_counts_true)
  pi_0_true <- sum(stayer_vec_true)/n

  all_seperated <- Seperator(data_obs)


  hpv_obs <- IntroduceMissing(all_seperated[[1]],.055)
  cyt_obs <- IntroduceMissing(all_seperated[[2]], .003)
  colpo_obs <- IntroduceMissingColpo(hpv_obs, cyt_obs, all_seperated[[3]])
  hpv_obs[hpv_obs > 1] <- 1
  hpv_pers_obs <- IntroducePersistence(hpv_obs)

  data_obs_list <- IntroduceMissingAll(hpv_pers_obs,cyt_obs, colpo_obs, c(0.001,0.157,0.206,0.228,0.169))
  return(data_obs_list)
}

#' Standard Initial vector
#'
#' Returns standard initial vector
#'
#' @return Initial Vector
#'
#' @export
GetInit <- function(){
  init_hc <- c(0.272,0.004,0.112,0.070,0.003,0.001,0.047,0.025,0.106,0.057,0.014,0.020,0.047,0.025,0.106,0.057,0.014,0.020)
  return(init_hc)
}

#' Standard Transition Matrix
#'
#' Returns standard transition matrix
#'
#' @return Transition matrix
#'
#' @export
GetTran <- function(){
  tran_hc <- matrix(c(0.809,0.000,0.089,0.000,0.002,0.000,0.073,0.000,0.026,0.000,0.001,0.0,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.000,0.130,0.000,0.309,0.000,0.000,0.000,0.426,0.000,0.136,0.000,0.0,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.643,0.000,0.260,0.000,0.007,0.000,0.043,0.000,0.045,0.000,0.002,0.0,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.000,0.021,0.000,0.730,0.000,0.002,0.000,0.128,0.000,0.119,0.000,0.0,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.554,0.000,0.290,0.000,0.071,0.000,0.031,0.000,0.054,0.000,0.000,0.0,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.000,0.000,0.000,0.500,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.5,0.000,0.000,0.000,0.000,0.000,0.000,
                      0.596,0.000,0.056,0.000,0.001,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.206,0.000,0.131,0.000,0.010,0.000,
                      0.000,0.043,0.000,0.121,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.509,0.000,0.304,0.000,0.023,
                      0.367,0.000,0.103,0.000,0.011,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.167,0.000,0.321,0.000,0.031,0.000,
                      0.000,0.015,0.000,0.145,0.000,0.002,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.274,0.000,0.499,0.000,0.065,
                      0.357,0.000,0.152,0.000,0.015,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.134,0.000,0.231,0.000,0.111,0.000,
                      0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.316,0.000,0.474,0.000,0.211,
                      0.435,0.000,0.043,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.320,0.003,0.186,0.002,0.011,0.000,
                      0.000,0.024,0.000,0.065,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.560,0.000,0.325,0.000,0.025,
                      0.285,0.000,0.080,0.000,0.005,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.187,0.002,0.395,0.003,0.042,0.000,
                      0.000,0.016,0.000,0.084,0.000,0.003,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.246,0.000,0.589,0.000,0.063,
                      0.207,0.000,0.136,0.000,0.007,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.121,0.001,0.405,0.003,0.118,0.001,
                      0.000,0.000,0.000,0.042,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.000,0.042,0.000,0.583,0.000,0.333),18,18, byrow = T)
  return(tran_hc)
}

#' Standard Classification Matrix
#'
#' Returns standard classification matrix
#'
#' @return Classification matrix
#'
#' @export
GetClass <- function(){
  class_hc <- matrix(0,18,18)
  for (i in 1:18){
    if ((i-1) %% 2 == 0){
      class_hc[i,i] <- .95
      if ((((i-1) %% 6) %/% 2) == 0){
        class_hc[i,i+2] <- (1 - class_hc[i,i])/2
        class_hc[i,i+4] <- (1 - class_hc[i,i])/2
      }  else if ((((i-1) %% 6) %/% 2) == 1){
        class_hc[i,i-2] <- (1 - class_hc[i,i])/2
        class_hc[i,i+2] <- (1 - class_hc[i,i])/2
      } else if ((((i-1) %% 6) %/% 2) == 2){
        class_hc[i,i-2] <- (1 - class_hc[i,i])/2
        class_hc[i,i-4] <- (1 - class_hc[i,i])/2
      }
    }else {
      class_hc[i,i] <- .95
      if ((((i-1) %% 6) %/% 2) == 0){
        class_hc[i,i-1] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i+1] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i+2] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i+3] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i+4] <- (1 - class_hc[i,i]) / 5
      } else if ((((i-1) %% 6) %/% 2) == 1){
        class_hc[i,i-3] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i-2] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i-1] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i+1] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i+2] <- (1 - class_hc[i,i]) / 5
      } else if ((((i-1) %% 6) %/% 2) == 2){
        class_hc[i,i-5] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i-4] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i-3] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i-2] <- (1 - class_hc[i,i]) / 5
        class_hc[i,i-1] <- (1 - class_hc[i,i]) / 5
      }
    }
  }
  return(class_hc)
}

##################################################
################ Method Functions ################
##################################################

#' Calculates State Given Individual Test Results
#'
#' Calculates state tuple result given individual tests, can output partially observed data. Incorportated in CombineandPattern().
#'
#' @param val1 value from test 1
#' @param val2 value from test 2
#' @param val3 value from test 3
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#' @param reduce_ng if true, if the first two tests are negative imputes negative for the third test
#'
#' @return tuple value, may be vector if partially observed
#'
#' @export
StateCalculator <- function(val1, val2, val3, max2, max3, reduce_ng = F){
  max1 <- 3
  to_return <- c()
  if (!is.na(val1) & !is.na(val2) & !is.na(val3)){
    return(((max2*max3) * val1) + (max3 * val2) + val3)
  } else if (!is.na(val1) & !is.na(val2) & is.na(val3)){
    if (reduce_ng){
      if(val1 == 0 & val2 == 0){
        return(0)
      }
    }
    for (i in 0:(max3-1)){
      to_return <- c(to_return,(((max2*max3) * val1) + (max3 * val2) + i))
    }
  } else if (!is.na(val1) & is.na(val2) & !is.na(val3)){
    for (i in 0:(max2 -1)){
      to_return <- c(to_return,(((max2*max3) * val1) + (max3 * i) + val3))
    }
  } else if (is.na(val1) & !is.na(val2) & !is.na(val3)){
    for (i in 0:(max1 -1)){
      to_return <- c(to_return,(((max2*max3) * i) + (max3 * val2) + val3))
    }
  } else if (!is.na(val1) & is.na(val2) & is.na(val3)){
    for (i in 0:(max2 -1)){
      for (j in 0:(max3 -1)){
        to_return <- c(to_return,(((max2*max3) * val1) + (max3 * i) + j))
      }
    }
  } else if (is.na(val1) & !is.na(val2) & is.na(val3)){
    for (i in 0:(max1 -1)){
      for (j in 0:(max3 -1)){
        to_return <- c(to_return,(((max2*max3) * i) + (max3 * val2) + j))
      }
    }
  } else if (is.na(val1) & is.na(val2) & !is.na(val3)){
    for (i in 0:(max1 -1)){
      for (j in 0:(max2 -1)){
        to_return <- c(to_return,(((max2*max3) * i) + (max3 * j) + val3))
      }
    }
  } else {return(NA)}
  return(to_return)
}

#' Combines three data matrices
#'
#' Combines three data matrices to one array
#'
#' @param mat1 first source data.  Must be array with a third dimension of at least 2.  3rd dimension can be all NAs.  Incorporated in CombineandPattern
#' @param mat2 second data matrix
#' @param mat3 third data matrix
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#' @param reduce T/F that is used as reduce_ng in StateCalculator
#'
#' @return combined array
#'
#' @export
Combiner <- function(mat1, mat2, mat3, max2, max3,reduce){
  max_part <- max(2 * max2 * max3)
  full_data <- array(NA, dim = c(dim(mat1)[1], dim(mat1)[2], max_part))
  for (i in 1:dim(mat1)[1]){
    for (j in 1:dim(mat1)[2]){
      vals <- c()
      num_of_parts <- sum(!is.na(mat1[i,j,]))

      if (num_of_parts == 0 | num_of_parts == 1){
        vals <- StateCalculator(mat1[i,j,1], mat2[i,j], mat3[i,j], max2, max3,reduce)
      } else {
        for (val_index in 1:num_of_parts){
          vals <- c(vals, StateCalculator(mat1[i,j,val_index], mat2[i,j], mat3[i,j], max2, max3,reduce))
          vals <- sort(unique(vals))
        }
      }

      for (part_val in 1:length(vals)){
        full_data[i,j,part_val] <- vals[part_val]
      }
    }
  }
  return(full_data)
}

#' Normalized Matrices
#'
#' Normalizes matrices so rows sum to one
#'
#' @param data matrix that will be row normalied
#'
#' @return row normalized matrix
#'
#' @export
Normalize <- function(data){
  for (i in 1:dim(data)[1]){
    if (sum(data[i,],na.rm = TRUE) != 0){
      data[i,] <- data[i,]/sum(data[i,])
    }
  }
  return(data)
}

#' Classifies observed given truth
#'
#' Returns classification probability of x_val as y_val given classification probability
#'
#' @param x_val observed value
#' @param y_val true value
#' @param class classification matrix
#'
#' @export
Classification <- function(x_val,y_val,class){
  if (is.na(x_val[1])){
    return(1)
  } else {
    rs <- 0
    for (val in x_val){
      if (is.na(val)){
        return(rs)
      } else {
        rs <- rs + class[[y_val + 1,val + 1]]
      }
    }
  }
  return(rs)
}

#' Returns Number of states
#'
#' Returns Number of states
#'
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#'
#' @return vector from 0 to k where k is the max value of the states
#'
#' @export
numOfStates <- function(max2,max3){
  max1 <- 3
  max_val <- max1 * max2 * max3
  return(c(0:(max_val-1)))
}

#' Patternizes Data
#'
#' Gets all patterns in data. Incorporated in CombineandPattern
#'
#' @param states large data array
#'
#' @return dataframe of patterns and frequencies.
#'
#' @export
GetPatterns <- function(states){
  options(scipen = 999)
  vals <- rep(NaN,dim(states)[1])
  for (i in 1:dim(states)[1]){
    val <- ""
    for(j in 1:dim(states)[2]){
      for(k in 1:dim(states)[3]){
        if (!is.na(states[i,j,k])){
          states_val <- states[i,j,k]
        } else {
          states_val <- NA
        }
        val <- paste0(val, states_val,",")
      }
    }
    val <- substr(val, 1, nchar(val)-1)
    vals[i] <- val
  }
  return(vals)
}

#' Turns patterns into data array
#'
#' Turns patterns into data array.  Incorporated in CombineandPattern
#'
#' @param unique_patterns patterns output from GetPatterns
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#' @param time_length number of observations per person
#'
#' @return list of array of unique patterns and frequency vector
#'
#' @export
Pattern2Data <- function(unique_patterns,max2,max3,time_length){
  max_part <- max(2 * max2 * max3)
  n <- length(unique_patterns[[1]])
  freq_vec <- numeric(n)
  state_array <- array(0, dim = c(n,time_length,max_part))
  for (i in 1:n){
    freq_vec[i] <- unique_patterns$Freq[[i]]
    pattern <- unique_patterns$all_patterns[[i]]

    pattern <- strsplit(pattern, split = ",")[[1]]
    for (time in 1:time_length){
      state_array[i,time,] <- suppressWarnings(as.integer(pattern[(((time-1)*max_part)+1):(time*max_part)]))
    }
  }

  return(list(state_array, freq_vec))
}

#' Adapted Forward Algorithm
#'
#' Adapted forward that incorporates partially observed data.  Partial data is resolved in linear time.
#'
#' @param data two dimensional matrix for individual.  1st dimension is time 2nd is partially observed data
#' @param time what time to caluclate for.  Works recursivly, if time = t calculates for all values less than t
#' @param init initial probabilities
#' @param tran transition probabilities
#' @param class classification probabilities
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#'
#' @return list (each individual) of list (each time) of list (each partially observed data) of vector (forward probabilities)
#'
#' @export
ForwardLinear <- function(data,time,init,tran,class,max2,max3){
  k <- numOfStates(max2, max3)

  alpha_matrix<-vector("list",time)
  alpha_i <- numeric(length(k))
  num_of_part_vec <- numeric(dim(data)[1])

  for (i in 1:time){
    num_of_part_vec[i] <- sum(!is.na(data[i,]))
    if (num_of_part_vec[i] == 0){num_of_part_vec[i] = 1}
    alpha_matrix[[i]] <- vector("list",num_of_part_vec[i])
    for (j in 1:num_of_part_vec[i]){
      alpha_matrix[[i]][[j]] <- alpha_i
    }
  }

  partial_vals <- data[1,][!is.na(data[1,])]
  for (i in 1:num_of_part_vec[1]){
    for (j in 1:length(k)){
      alpha_matrix[[1]][[i]][j] <- init[j] * Classification(partial_vals[i],k[j],class)
    }
  }

  if (time > 1){
    for (i in 2:time){
      partial_vals <- data[i,][!is.na(data[i,])]

      if (num_of_part_vec[i - 1] == 1){
        old_alpha <- alpha_matrix[[i-1]][[1]]
      } else {
        old_alpha <- numeric(length(k))
        for (j in 1:num_of_part_vec[i-1]){
          old_alpha <- old_alpha + alpha_matrix[[i-1]][[j]]
        }
      }

      for (j in 1:num_of_part_vec[i]){
        for (l in 1:length(k)){
          alpha_matrix[[i]][[j]][l] <- (old_alpha %*% tran[,l]) *  Classification(partial_vals[j],k[l],class)
        }
      }
    }
  }

  return(alpha_matrix)
}

#' Adapted Backward Algorithm
#'
#' Adapted backward that incorporates partially observed data.  Partial data is resolved in linear time.
#'
#' @param data two dimensional matrix for individual.  1st dimension is time 2nd is partially observed data
#' @param time what time to caluclate for.  Works recursivly, if time = t calculates for all values greater than t
#' @param tran transition probabilities
#' @param class classification probabilities
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#'
#' @return list (each individual) of list (each time) of list (each partially observed data) of vector (forward probabilities)
#'
#' @export
BackwardLinear <- function(data,time,tran,class,max2,max3){
  k <- numOfStates(max2,max3)

  beta_i <- numeric(length(k))
  misclass_vector <- numeric(length(k))
  beta_matrix<-vector("list",dim(data)[1])
  num_of_part_vec <- numeric(dim(data)[1])


  for (i in 1:dim(data)[1]){
    num_of_part_vec[i] <- sum(!is.na(data[i,]))
    if (num_of_part_vec[i] == 0){num_of_part_vec[i] = 1}
  }


  for (i in 1:(dim(data)[1] - 1)){
    beta_matrix[[i]] <- vector("list",num_of_part_vec[i + 1])
    for (j in 1:num_of_part_vec[i + 1]){
      beta_matrix[[i]][[j]] <- beta_i
    }
  }

  beta_matrix[[dim(data)[1]]] <- vector("list",1)
  for (i in 1:length(k)){
    beta_matrix[[dim(data)[1]]][[1]][i] <- 1
  }

  if (time != dim(data)[1]) {
    for (i in (dim(data)[1]-1):time){
      partial_vals <- data[i+1,][!is.na(data[i+1,])]

      new_beta <- numeric(length(k))
      for (m in 1:length(beta_matrix[[i+1]])){
        new_beta <- new_beta + beta_matrix[[i+1]][[m]]
      }


      for (m in 1:num_of_part_vec[i+1]){

        for (l in 1:length(k)){
          misclass_vector[l] <- Classification(partial_vals[m],k[l],class)
        }

        for (j in 1:length(k)){
          beta_matrix[[i]][[m]][j] <- sum(new_beta * tran[j,] * misclass_vector)
        }

      }
    }
  }


  return(beta_matrix)
}

#' Calculates Probability Latent Values are 0 given observed data
#'
#' Calculates probability for individuals observed values assuming latent values are zero
#'
#' @param data dimensional matrix for individual.  1st dimension is time 2nd is partially observed data
#' @param class classification matrix
#'
#' @return probability of obsevred values if latent values are all 0
#'
#' @export
ProdClass <- function(data,class){
  prod <- 1
  for (i in 1:dim(data)[1]){
    prod <- prod * Classification(data[i,],0,class)
  }
  return(prod)
}

#' Calculates Likelihood of Movers
#'
#' Calculates likelihood of movers by summing over forward quantity
#'
#' @param forw forward quantity for every individual
#' @param pi probability of being a stayer
#'
#' @return mover likelihood
#'
#' @export
CalcLikelihoodMover <- function(forw, pi){

  rs <- 0
  for (i in 1:length(forw[[length(forw)]])){
    rs <- rs + forw[[length(forw)]][[i]]
  }

  return(sum(rs)[[1]]* (1 - pi))
}

#' Calculates Likelihood of Stayers
#'
#' Calculates Likelihood of Stayers
#'
#' @param data two dimensional matrix for individual.  1st dimension is time 2nd is partially observed data
#' @param class classification matrix
#' @param pi probability of being a stayer
#'
#' @return stayer likelihood
#'
#' @export
CalcLikelihoodStayer <- function(data, class, pi){
  return((pi * ProdClass(data,class)))
}

#' Estimates Stayer Proportion
#'
#' Calculates first expectation in supplementary materials section 1.2
#'
#' @param data two dimensional matrix for individual.  1st dimension is time 2nd is partially observed data
#' @param init initial probabilities
#' @param tran transition probabilities
#' @param class classification probabilities
#' @param freq_vec frequency vector to go along with data
#' @param pi_0 probability of being stayer
#' @param likelihoods list of likelihoods for each individual
#'
#' @return estimates stayer probability
#'
#' @export
CalcStayerLin <- function(data,init,tran,class,freq_vec,pi_0, likelihoods){

  CalcStayerHelper <- function(data,likelihood){
    return((pi_0 * ProdClass(data,class) / likelihood))
  }


  pi_vec <- numeric(dim(data)[1])
  for (i in 1:dim(data)[1]){
    pi_vec[i] <- CalcStayerHelper(data[i,,], likelihoods[i])
  }
  pi_0 <- (pi_vec %*% freq_vec)/sum(freq_vec)
  return(pi_0[1])
}

#' Estimates Initial Probabilities
#'
#' Calculates second expectation in supplementary materials section 1.2
#'
#' @param data two dimensional matrix for individual.  1st dimension is time 2nd is partially observed data
#' @param freq_vec frequency vector to go along with data
#' @param pi_0 probability of being stayer
#' @param forw forward quantity for each individual
#' @param backw backward quantity for each individual
#' @param likelihoods list of likelihoods for each individual
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#'
#' @return estimates vector of initial probabilities
#'
#' @export
CalcInitialLin <- function(data,freq_vec, pi_0, forw, backw, likelihoods,max2,max3){
  k <- numOfStates(max2,max3)
  prob_list <- numeric(length(k))

  for (i in 1:dim(data)[1]){

    ind_likelihood <- likelihoods[[i]]
    forw_sum <- numeric(length(k))
    backw_sum <- numeric(length(k))

    for (j in 1:length(forw[[i]][[1]])){
      forw_sum <- forw_sum + forw[[i]][[1]][[j]]
    }

    for (j in 1:length(backw[[i]][[1]])){
      backw_sum <- backw_sum + backw[[i]][[1]][[j]]
    }

    mover <- forw_sum * backw_sum * (1 - pi_0) / ind_likelihood
    prob_list <- prob_list + (mover * freq_vec[i])
  }
  return(prob_list/sum(prob_list))
}


#' Estimates Transition Probabilities
#'
#' Calculates third expectation in supplementary materials section 1.2
#'
#' @param data two dimensional matrix for individual.  1st dimension is time 2nd is partially observed data
#' @param tran transition probabilities
#' @param class classification probabilities
#' @param freq_vec frequency vector to go along with data
#' @param pi_0 probability of being stayer
#' @param forw forward quantity for each individual
#' @param backw backward quantity for each individual
#' @param likelihoods list of likelihoods for each individual
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#'
#' @return estimates matrix of transition probabilities
#'
#' @export
CalcTransitionLin <- function(data, tran,class,freq_vec, pi_0, forw, backw, likelihoods,max2,max3){

  k <- numOfStates(max2,max3)
  tran_matrix <- matrix(0L, nrow = length(k), ncol = length(k))

  for (ind in 1:dim(data)[1]){
    forward <- forw[[ind]]
    backward <- backw[[ind]]
    denom <- likelihoods[[ind]]
    for (time in 2:dim(data)[2]){


      forw_sum <- numeric(length(k))
      backw_sum <- numeric(length(k))

      for (i in 1:length(forward[[time-1]])){
        forw_sum <- forw_sum + forward[[time-1]][[i]]
      }

      for (i in 1:length(backward[[time]])){
        backw_sum <- backw_sum + backward[[time]][[i]]
      }


      tran_matrix_current <- matrix(0L, nrow = length(k), ncol = length(k))

      for (initial_state in 1:length(k)){
        for(new_state in 1:length(k)){
          num <- forw_sum[[initial_state]] * backw_sum[[new_state]] * tran[initial_state,new_state] * Classification(data[ind,time,],new_state-1,class) * (1 - pi_0)
          if(denom != 0){
            tran_matrix_current[initial_state,new_state] <- tran_matrix_current[initial_state,new_state] + (num/denom)
          }
        }
      }
      tran_matrix <- tran_matrix + ((tran_matrix_current) * freq_vec[ind])
    }
  }
  return(Normalize(tran_matrix))
}

#' Estimates Classification Probabilities
#'
#' Calculates fourth expectation in supplementary materials section 1.2
#'
#' @param data two dimensional matrix for individual.  1st dimension is time 2nd is partially observed data
#' @param class classification probabilities
#' @param freq_vec frequency vector to go along with data
#' @param pi_0 probability of being stayer
#' @param forw forward quantity for each individual
#' @param backw backward quantity for each individual
#' @param likelihoods list of likelihoods for each individual
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#'
#' @return estimates matrix of classification probabilities
#'
#' @export
CalcClassificationLin <- function(data,class, freq_vec, pi_0, forw, backw, likelihoods,max2,max3){
  k <- numOfStates(max2,max3)
  num_matrix_full <- matrix(0L, nrow = length(k), ncol = length(k))
  for (ind in 1:dim(data)[1]){

    forward <- forw[[ind]]
    backward <- backw[[ind]]
    individual_likelihood <- likelihoods[[ind]]

    for(time in 1:dim(data)[2]){
      num_matrix <- matrix(0L, nrow = length(k), ncol = length(k))
      backward_sum <- numeric(length(k))

      for (i in 1:length(backward[[time]])){
        backward_sum <- backward_sum + backward[[time]][[i]]
      }


      x_vals <- data[ind,time,][!is.na(data[ind,time,])]

      for (u in 1:length(forward[[time]])){
        x_val <- x_vals[u]
        if (!is.na(x_val)){


          chain <- data[ind,,]
          chain[time,] <- NA
          chain[time,1] <- x_val

          mover <- (forward[[time]][[u]]*backward_sum * (1 - pi_0) / individual_likelihood)

          stayer <- (pi_0 * ProdClass(chain,class)) / individual_likelihood

          mover_stayer <- mover
          mover_stayer[1] <- mover_stayer[1] + stayer

          for (y_state in 1:length(k)){
            num_matrix[y_state,x_val + 1] <- num_matrix[y_state, x_val+1] + mover_stayer[y_state]
          }
        }
      }


      #Normalization for partial
      if (sum(num_matrix) != 0){
        num_matrix <- num_matrix/sum(num_matrix)
      }
      num_matrix_full <- num_matrix_full + (num_matrix * freq_vec[ind])

    }
  }

  return(Normalize(num_matrix_full))
}

#' Runs EM Algorithm
#'
#' Runs EM using adapted forward-backward that calculates partially observed data.
#'    Calculates probability of being a stayer, intiial vector, transition matrix, and classification matrix
#'    Uses the following functions: ForwardLinear, BackwardLinear, CalcLikelihoodMover, CalcLikelihoodStayer,
#'        CalcInitialLin, CalcTransitionLin, CalcClassificationLin, and CalcStayerLin
#'
#' @param data_array data array of three combined tests, patternized
#' @param freq_vec frequency vector of patterns
#' @param epsilon threshold for EM convergence, stops when likelihood percentage increase is below epsilon
#' @param init initial probabilities
#' @param tran transition probabilities
#' @param class classification probabilities
#' @param pi_0 probability of being stayer
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#' @param time_length Number of observations per individual
#'
#' @returns list of initial parameters, estimated parameters and likelihood
#'     initial and estimated parameter entries have same format
#'     1) initial 2) transition 3) classification 4) stayer probability
#'
#' @export
EM <- function(data_array, freq_vec, epsilon, init, tran, class, pi_0, max2, max3, time_length){

  init_og <- init
  tran_og <- tran
  class_og <- class
  pi_0_og <- pi_0

  forw <- apply(data_array,1,ForwardLinear, time = time_length, init,tran,class,max2,max3)
  backw <- apply(data_array,1,FUN = BackwardLinear, time = 1, tran, class,max2,max3)
  likelihoodsMover <- unlist(lapply(forw, CalcLikelihoodMover, pi = pi_0))
  likelihoodsStayers <- apply(data_array,1,FUN = CalcLikelihoodStayer, class, pi_0)
  likelihoods <- likelihoodsStayers + likelihoodsMover
  old_likelihood <- (log(likelihoods) %*% freq_vec)[1]

  init_new <- CalcInitialLin(data_array,freq_vec,pi_0,forw,backw,likelihoods,max2,max3)
  tran_new <- CalcTransitionLin(data_array,tran,class,freq_vec,pi_0,forw,backw,likelihoods,max2,max3)
  class_new <- CalcClassificationLin(data_array,class,freq_vec,pi_0,forw,backw,likelihoods,max2,max3)
  pi_0_new <- CalcStayerLin(data_array,init,tran,class,freq_vec,pi_0,likelihoods)

  init <- init_new
  tran <- tran_new
  class <- class_new
  pi_0 <- pi_0_new

  forw <- apply(data_array,1,ForwardLinear, time = time_length, init,tran,class,max2,max3)
  backw <- apply(data_array,1,FUN = BackwardLinear, time = 1, tran, class,max2,max3)
  likelihoodsMover <- unlist(lapply(forw, CalcLikelihoodMover, pi = pi_0))
  likelihoodsStayers <- apply(data_array,1,FUN = CalcLikelihoodStayer, class, pi_0)
  likelihoods <- likelihoodsStayers + likelihoodsMover
  new_likelihood <- (log(likelihoods) %*% freq_vec)[1]

  print(-(new_likelihood - old_likelihood)/old_likelihood)

  while (-(new_likelihood - old_likelihood)/old_likelihood > epsilon){

    old_likelihood <- new_likelihood

    init_new <- CalcInitialLin(data_array,freq_vec,pi_0,forw,backw,likelihoods,max2,max3)
    tran_new <- CalcTransitionLin(data_array,tran,class,freq_vec,pi_0,forw,backw,likelihoods,max2,max3)
    class_new <- CalcClassificationLin(data_array,class,freq_vec,pi_0,forw,backw,likelihoods,max2,max3)
    pi_0_new <- CalcStayerLin(data_array,init,tran,class,freq_vec,pi_0,likelihoods)

    init <- init_new
    tran <- tran_new
    class <- class_new
    pi_0 <- pi_0_new

    forw <- apply(data_array,1,ForwardLinear, time = time_length, init,tran,class,max2,max3)
    backw <- apply(data_array,1,FUN = BackwardLinear, time = 1, tran, class,max2,max3)
    likelihoodsMover <- unlist(lapply(forw, CalcLikelihoodMover, pi = pi_0))
    likelihoodsStayers <- apply(data_array,1,FUN = CalcLikelihoodStayer, class, pi_0)
    likelihoods <- likelihoodsStayers + likelihoodsMover
    new_likelihood <- (log(likelihoods) %*% freq_vec)[1]

    print(-(new_likelihood - old_likelihood)/old_likelihood)
  }
  initial_parameters <- list(init_og,tran_og, class_og, pi_0_og)
  estimated_parameters <- list(init, tran, class, pi_0)
  to_save <- list(initial_parameters,estimated_parameters, new_likelihood)
  return(to_save)
}

#' Combines three tests and patternizes them
#'
#' Combines three tests and patternizes them.  Matrix 2 and 3 must be of the same size.  Matrix 1 must be an array with
#'     the first two dimensions matching those of matrix 2 and 3 and with a third dimension of atleast 2.
#'
#' @param mat1_pers 3 dimensional data array of test 1
#' @param mat2 matrix for test 2
#' @param mat3 matrix for test 3
#' @param max2 number of values the second test can take on
#' @param max3 number of values the third test can take on
#' @param reduce for StateCalculator if true, if the first two tests are negative imputes negative for the third test in Combiner
#' @param time_length Number of observations per individual
#'
#' @return (list of patternized data and frequency vector)
#'
#' @export
CombineandPattern <- function(mat1_pers,mat2,mat3, max2, max3, reduce,time_length){
  combined_data <- Combiner(mat1_pers,mat2,mat3, max2, max3,reduce)
  all_patterns <- GetPatterns(combined_data)
  unique_patterns <- as.data.frame(table(all_patterns), stringsAsFactors = F)
  data_pattern_list <- Pattern2Data(unique_patterns,max2,max3,time_length)

  data_pattern <- data_pattern_list[[1]]
  freq_vec <- data_pattern_list[[2]]

  return(data_pattern_list)
}


############

# max2 <- 3
# max3 <- 2
#
# three_mats <- GenerateSimulatedData(3500,5,.05,max2,max3)
#
# init <- GetInit()
# tran <- GetTran()
# class <- GetClass()
# pi_0 <- .05
#
# data_pattern_list <- CombineandPattern(three_mats[[1]],three_mats[[2]],three_mats[[3]],max2,max3,T,5)
# data_pattern <- data_pattern_list[[1]]
# freq_vec <- data_pattern_list[[2]]
#
# parameters <- EM(data_pattern,freq_vec,.01,init,tran,class,pi_0,max2,max3,5)
