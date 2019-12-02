## R package 'PartiallyObservedHMM'

Creates hidden Markov model (HMM) that accounts for partially observed data.  Partially observed data is represented as a 3-dimensional array where all possibilities for individual "i" at time "t" can be visualized as a vector in row "i" column "t" of the 3-d array.  For our model, data must be generated from three unique tests.  The first data source must be a trichotomous 3-d array (that may or may not include partially observed data).  The second and third data sources must be 2-dimensional matrices. The second and third data sources must be identical in dimension along with the first two dimensions of the first data source.  

In our example the first data source is an HPV test.  Persistence has been added so the data is now trichotomous (negative, newly-positive, persistent).  The second source is a pap smear, and the final source is a colposcopy.  

-----

## 1\. Installation 

```{r}
install.packages("devtools")
devtools::install_github("jordanaron22/PartiallyObservedHMM")
```

-----
