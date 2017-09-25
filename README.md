#README
-----------

This package is related to the work in "Sparse Regularization in Marketing and Economics" by Guanhao Feng, Nicholas Polson, Yuexi Wang and Jianeng Xu. <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3022856>

##Functions

The package contains five functions and their usage is similar to those in `glmnet` package.

- alphanorm: build the alpha-norm model with input data. The coefficients are solved via a proximal algorithm and coordinate descent.

- coef.alphanorm: take an alphanorm object as input and output the coefficients in the object

- cv.alphanorm: use cross-validation to find the best tuning parameter $\alpha$ and $\lambda$ 

- plot.alphanorm: plot the coefficient profile with respect to either $log(\lambda)$, or $l_\alpha$ norm of the coefficients 

- predict.alphanorm: get predicted value from a new input and a fitted alphanorm model
   
Detailed usage of functions can be found in [alphanorm.pdf](https://github.com/yxwang99/alphanorm/blob/master/alphanorm.pdf)  
   
## Installation of Package

The package is currently not available on R CRAN. You can use it via:

```
devtools::install_github(yxwang99/alphanorm)
library(alphanorm)
```

If you encounter any problem using this package, please email to yxwang99@uchicago.edu





