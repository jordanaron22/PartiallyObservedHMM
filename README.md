## R package `PartiallyObservedHMM`

This package creates a hidden Markov model (HMM) that accounts for partially observed data.  Partially observed data is represented as a 3-dimensional array where all possibilities for individual "i" at time "t" can be visualized as a vector in row "i" column "t" of the 3-d array.  For our model, data must be generated from three unique tests.  The first data source must be a trichotomous 3-d array (that may or may not include partially observed data).  The second and third data sources must be 2-dimensional matrices. The second and third data sources must be identical in dimension along with the first two dimensions of the first data source.  

In our example the first data source is an HPV test.  Persistence has been added so the data is now trichotomous (negative, newly-positive, persistent).  The second source is a pap smear, and the final source is a colposcopy.  

-----
## News

- 12/02/19: Added 'read.me'
- 02/10/20: Added Rafferty method functions and documentation 

-----

## 1\. Installation 

```{r}
install.packages("devtools")
devtools::install_github("jordanaron22/PartiallyObservedHMM")
```
-----

## 2\. Examples 

### 2.1 First Order 

#### 2.1.1 Simulated Data 

```{r}
three_mats <- GenerateSimulatedData(n,t,pi_0)

data_pattern_list <- CombineandPattern(three_mats[[1]],three_mats[[2]],three_mats[[3]],reduce,t)
data_pattern <- data_pattern_list[[1]]
freq_vec <- data_pattern_list[[2]]

init <- GetInit()
tran <- GetTran()
class <- GetClass()
pi_0 <- .05
max2 <- 3
max3 <- 2

parameters <- EM(data_pattern,freq_vec,epsilon,t,init,tran,class,pi_0,max2,max3)
```

Where n is the number of individuals, t is the number of observations per individual, p is the proportion of stayers, reduce is a T/F for if the third test value should be imputed to 0 if the first two test values are 0.

init, tran, class, and pi_0 are the parameters that we are optimizing.  Let k be the number of states.  init is a k length vector giving the probability of initially starting in any state.  tran is a k x k matrix denoting the probability of transitioning from one state to another (from row to column).  class is a k x k matrix with the probabilities of classifying the latent (true) state as the observed state (latent is row, observed is column).  pi_0 is the probability of being a stayer.  max2 is the number of outcomes of the second test, and max3 is for the third test.

epsilon is the likelihood threshold for EM convergence. 

#### 2.1.2 Other Data 

```{r}
data_pattern_list <- CombineandPattern(mat_one,mat_two,mat_three,reduce,t, max2, max3)
data_pattern <- data_pattern_list[[1]]
freq_vec <- data_pattern_list[[2]]

parameters <- EM(data_pattern, freq_vec, epsilon, t, initial, transition, classification, pi_0, max2, max3)
```

Where the variables are defined as above and: mat_one, mat_two, and mat_three are data sources for three unique tests, max2 is the number of possible outcomes in mat_two, max2 is the number of possible outcomes in mat_three, initial is an initial estimate for the vector of initial state probabilities, transition is an initial estimate for the state transition matrix, classification is an initial estimate for the state classification matrix, and pi_0 is an initial estimate of the proportion of stayers. 

### 2.2 Second Order (Rafferty Method)

#### 2.2.1 Simulated Data

```{r}
three_mats <- GenerateSimulatedDataRaff(n, t, pi_0, lambda)
data_pattern_list <- CombineandPattern(three_mats[[1]][[1]], three_mats[[1]][[2]], three_mats[[1]][[3]], reduce, t)
data_pattern <- data_pattern_list[[1]]
freq_vec <- data_pattern_list[[2]]

init <- GetInitRaff()
tran <- GetTran()
class <- GetClass()
pi_0 <- .05
lambda <- c(.7,.3)
max2 <- 3
max3 <- 2

parameters <- EMRaff(data_pattern, freq_vec, epsilon, t, init, tran, class, pi_0, lambda, max2, max3)
```

Variables are defined similarly to 2.1.1, except for init and lambda.  init is now a k by k matrix (where k is the number of states) giving the probability of being in the first two states (the row determines the first state and the column determines the second).  lambda is a vector of length two that gives the weights for the first and second order transition from the Rafferty method (the first value in the vector corresponds to the first order weight.  

#### 2.2.2 Other Data

```{r}
data_pattern_list <- CombineandPattern(mat_one, mat_two, mat_three, reduce, t, max2, max3)
data_pattern <- data_pattern_list[[1]]
freq_vec <- data_pattern_list[[2]]

parameters <- EMRaff(data_pattern,freq_vec,epsilon,t, initial, transition, classification, pi_0, lambda, max2, max3)
```

This is extremeley similar to 2.1.2, except for the use of EMRaff and the inclusion of the lambda variable
