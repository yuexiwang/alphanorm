# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#########
# code for alpha norm function
#########
#' @import graphics
#' @import stats

alpha_map = function(b,h,z,lambda,q,beta)
{
  tau = 0
  abz = abs(z)
  if(abz==h)
  {
    tau = ifelse(beta==0,0,sign(z)*b)
    return(tau)
  }

  if(abs(z)>h)
  {
    tau = sign(z)*uniroot(function(x) x+lambda*q*x^(q-1)-abz,
                          lower = b,
                          upper = abz)$root
    return(tau)
  }
  return(tau)
}

#' fit a sparse model with alpha-norm regularization
#' @description Fit a alph-norm model with proximal algorithm and coordinate descent
#' @param x the design matrix
#' @param y the response vector
#' @param lambda a vector of lambda values, default as exp(10:-10)
#' @param q a numerical value for q, 0<q<=1, with default 0.5
#' @param intercept whether the intercept term should be included, TRUE to be included(default), FALSE not to
#' @param tol tolerence of convergence condition
#' @param T number of maximum iterations for each coefficient
#' @param nlambda number of lambda wanted
#' @param trace print the process
#' @return An object of S3 class "alphanorm"
#' \item{\code{x}}{input design matrix}
#' \item{\code{y}}{input of response vector}
#' \item{\code{Lambda}}{input of lambda in the decreasing order}
#' \item{\code{q}}{input value of q}
#' \item{\code{Coefficient}}{matrix coefficients}
#' \item{\code{Intercept}}{non-penalized intercept(if intercept=TRUE), otherwise, NULL}
#' \item{\code{df}}{number of nonzero coefficients for each value of lambda}
#' @details
#' The sequence of models implied by \code{lambda} is solved via coordinate descent.
#' The objective function is:
#' \deqn{J(\beta)=1/2 \mbox{RSS}+\lambda*\mbox{penalty}}
#' Here the penalty is the \eqn{l_q} norm of coefficients, which is \eqn{\sum(|\beta_i|^q), 0<q<=1},
#' when \eqn{q=1}, it is actually same as lasso
#' @examples
#' x<-matrix(rnorm(100*100),100,100)
#' # Only the first 10 are true predictors
#' y<-x[,1:10]%*%rep(1,10)
#'
#' # Build a alpha-norm model
#' alphanorm.obj<-alphanorm(x,y,intercept=FALSE)
#' # Get coefficients
#' coef(alphanorm.obj)
#' # Get fitted values
#' predict(alphanorm.obj)
#' # Cross-validation to choose q and lambda
#' cv.alphanorm(x,y,intercept=FALSE)
#' # Plot coefficient profile according to log-lambda
#' plot(alphanorm.obj)
#' @references
#' Feng, Guanhao and Polson, Nick and Wang, Yuexi and Xu, Jianeng,
#' Sparse Regularization in Marketing and Economics (August 20, 2017).
#' Available at SSRN: \url{https://ssrn.com/abstract=3022856}
#'
#' Marjanovic, G. and V. Solo (2014). lq sparsity penalized linear regression with cyclic descent.
#' IEEE Transactions on Signal Processing 62(6), 1464–1475.
#' @seealso \code{\link{predict.alphanorm}}, \code{\link{coef.alphanorm}}, \code{\link{cv.alphanorm}}, and \code{\link{plot.alphanorm}} methods
#' @export
alphanorm = function(x,y,lambda=exp(10:-10),q=0.5,
                     intercept=TRUE,tol=1e-7,T=500,nlambda=NULL,trace=FALSE){


  if(!is.null(nlambda)){
    lambda<-exp(seq(10,-10,length.out = nlambda))
  }

  lambda = sort(lambda,decreasing = TRUE)
  L = length(lambda)
  N = dim(x)[1]
  P = dim(x)[2]
  y.mean = mean(y)
  x.mean = apply(x,2,mean)
  Beta.final = matrix(0,P,L)
  Inter = NULL

  # normalize
  if(intercept)
  {
    y = y - mean(y)
    x = sweep(x,2,x.mean)
  }
  x.norm = apply(x,2,function(x) sqrt(sum(x^2)))

  x = scale(x,center = FALSE,scale = x.norm)
  for(n_c in 1:dim(x)[2]){
    if(sum(is.na(x[,n_c]))>0){
      x[,n_c]<-0
    }
  }

  # initialize
  z.initial.max = max(abs(crossprod(x,y)))
  beta = rep(0,P)
  r = y

  # convergence condition
  Tol = tol*sum(r^2)/2
  Time = P*T
  b = (2*lambda*(1-q))^(1/(2-q))
  h = b+lambda*q*b^(q-1)

  for(l in 1:L)
  {
    if(trace==TRUE){
    cat("Computing: Log(lambda) = ",log(lambda[l]),"\n")
    }
    # proximal parameters and iteration initialization

    if(z.initial.max > h[l])
    {
      k = 1
      i = 1
      J.change.max = Tol*2
      J.change = rep(J.change.max,P)

      while(k < Time & J.change.max > Tol)
      {
        # coordinate
        if(i > P) i = 1

        # coordinate descent
        z = sum(x[,i]*r) + beta[i]
        tau = alpha_map(b[l],h[l],z,lambda[l],q,beta[i])
        r.change = (tau - beta[i])*x[,i]
        J.change[i] = sum((2*r - r.change)*r.change)/2 +
          lambda[l]*(abs(beta[i])^q-abs(tau)^q)

        # update beta, residual, objective function change
        r = r - r.change
        beta[i] = tau
        J.change.max = max(J.change)

        # next coordinate
        i = i+1
        k = k+1
      }


      if(k == Time & trace==TRUE)
        print("Coordinate descent did not converge. You might want to increase T or decrease tol")
      Beta.final[1:P,l] = beta
    }
  }

  Beta.final = Beta.final/x.norm
  Beta.final[which(!is.finite(Beta.final))]<-0
  if(intercept)
    Inter = y.mean - crossprod(x.mean,Beta.final)

  df=apply(Beta.final,2,function(x) sum(x!=0))
  alphanorm.obj<-list(x=x,y=y,"Lambda" = lambda,q=q,"Coefficient" = Beta.final,"Intercept" = Inter,df=df)
  class(alphanorm.obj)<-"alphanorm"

  return(alphanorm.obj)
}

#' Output the coefficients of "alphanorm" object
#' @param alphanorm.obj a fitted "alphanorm" object
#' @return coefficients of "alphanorm" object
#' @seealso \code{\link{alphanorm}}
#' @export
coef.alphanorm<-function(alphanorm.obj){
  if(is.null(alphanorm.obj$Intercept)){
    return(alphanorm.obj$Coefficient)
  }else{
    return(rbind(alphanorm.obj$Intercept,alphanorm.obj$Coefficient))
  }
}


#' Predict method for alpha-norm fits
#' @description Similar to other predict methods, this function predicts fitted values from a fitted alphanorm model
#' @param alphanorm.obj a fitted alpha-norm model, returned by alphanorm()
#' @param newx matrix of new values of x, if NULL, use the x in alphanorm.obj
#' @return  matrix of fitted values from alpha-norm model
#' @seealso \code{\link{alphanorm}}, and \code{\link{cv.alphanorm}} methods
#' @export
predict.alphanorm<-function(alphanorm.obj,newx=NULL){
  #The input is from the result of alpha norm

  if(is.null(newx)){
    newx<-alphanorm.obj$x
  }

  if(is.null(alphanorm.obj$Intercept)){
    yhat<-newx%*%alphanorm.obj$Coefficient
  }else{
    yhat<-cbind(1,newx)%*%rbind(alphanorm.obj$Intercept,alphanorm.obj$Coefficient)
  }

  return(yhat)
}

#' Cross-validation for alpha-norm
#' @description  Does k-fold cross-validation for alpha-norm, and return the best lambda and q
#' @param x design matrix
#' @param y response vector
#' @param lambda_Tune user-supplied lambda sequence
#' @param q_Tune user-supplied q sequence
#' @param intercept whether intercept should be in the model, default to be TRUE
#' @param nfolds number of folds , default to be 5
#' @param tol tolerence of convergence condition
#' @param T number of maximum iterations for each coefficient
#' @param trace print the process of alphanorm
#' @return An object of S3 class "cv.alphanorm"
#' \item{lambda}{the values of lambda used in the fits in the decreasing order}
#' \item{q}{the values of q used in the fits}
#' \item{cvm}{The mean cross-validation error, a matrix of length(q)*length(lambda) }
#' \item{lambda.min}{value of lambda that gives minimum cvm}
#' \item{q.min}{value of q that gives minimum cvm}
#' @seealso \code{\link{alphanorm}}
#' @export
cv.alphanorm<-function(x,y,lambda_Tune=exp(10:-10),q_Tune=c(0.1,0.5,0.9),
                       intercept=TRUE,nfolds=5,tol=1e-7,T=500,trace=FALSE){
  #Here we use mse as the measure for CV

  numTrain<-length(y)
  lambda_Tune<-sort(lambda_Tune,decreasing = TRUE)

  alphaNorm_mse<-array(NA,dim=c(length(q_Tune),length(lambda_Tune),nfolds))
  CVsample <- sample(1:numTrain,numTrain)
  CVcutoff <- round(seq(0, numTrain, length.out = nfolds + 1))

  for(j in 1:nfolds){
    newdata <- CVsample[(CVcutoff[j] + 1):CVcutoff[j + 1]]
    size <- length(newdata)
    X <- x[-newdata,]
    Y <- y[-newdata]
    Xnew <- x[newdata,]
    Ynew <- y[newdata]

    for(k in 1:length(q_Tune)){
      tmp_tune<-alphanorm(X,Y,lambda=lambda_Tune,q=q_Tune[k],intercept=intercept,trace=trace)
      tmp_pred<-predict(tmp_tune,newx=Xnew)
      tmp_mse<-apply(tmp_pred,2,function(x) mean((Ynew-x)^2))
      alphaNorm_mse[k,,j]<-tmp_mse
    }
  }

  alphaNorm_cve <- apply(alphaNorm_mse, c(1,2), mean)
  alphaNorm_best <- which(alphaNorm_cve == min(alphaNorm_cve),arr.ind=TRUE)

  best_lambda<- lambda_Tune[alphaNorm_best[2]]
  best_q <- q_Tune[alphaNorm_best[1]]

  cv.alphanorm.obj<-list(lambda=lambda_Tune,q=q_Tune,lambda.min=best_lambda,q.min=best_q,cvm=alphaNorm_cve)
  class(cv.alphanorm.obj)<-"cv.alphanorm"

  return(cv.alphanorm.obj)
}

#' plot coefficient for a "alphanorm"
#' @description Produce a coefficient profile plot of the coefficient paths for a fitted "alphanorm" object
#' @param alphanorm.obj fitted "alphanorm" model
#' @param xvar what is on the X-axis. "norm" plots against the $L_q$-norm of the coefficients, "lambda" against the log-lambda sequence
#' @param legend whether legend should be plotted
#' @examples
#' x=matrix(rnorm(100*20),100,20)
#' y=rnorm(100)
#' obj1=alphanorm(x,y)
#' plot(obj1)
#' plot(obj1,xvar="norm")
#' @seealso \code{\link{alphanorm}}
#' @export
plot.alphanorm<-function(alphanorm.obj,xvar=c("lambda"),legend=FALSE){

  if(xvar=="norm"){
    x<-apply(alphanorm.obj$Coefficient,2,function(x) sum(abs(x)^alphanorm.obj$q))
    matplot(x,t(alphanorm.obj$Coefficient),type="l",xlab="Lq-norm",ylab="Coefficients")
  }else if(xvar=="lambda"){
    x<-log(alphanorm.obj$Lambda)
    matplot(x,t(alphanorm.obj$Coefficient),type="l",xlab="log-lambda",ylab="Coefficients")
  }else{
    stop("X-axis can only be norm or lambda")
  }

  if(legend){
  if(xvar=="norm"){
    if(is.null(colnames(alphanorm.obj$x))){
    legend("topleft",legend=1:dim(alphanorm.obj$x)[2],col=1:6,lty=1:5)
    }else{
    legend("topleft",legend=colnames(alphanorm.obj$x),col=1:6,lty=1:5)
    }
  }else{
    if(is.null(colnames(alphanorm.obj$x))){
      legend("topright",legend=1:dim(alphanorm.obj$x)[2],col=1:6,lty=1:5)
    }else{
      legend("topright",legend=colnames(alphanorm.obj$x),col=1:6,lty=1:5)
    }
  }
  }
}



