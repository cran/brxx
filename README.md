---
title: "brxx"
output:
  html_document:
    toc: yes
pagetitle: brxx
---
Bayesian reliability estimation

Description
-----------
When samples contain missing data, are small, or are suspected of bias, estimation of scale reliability may not be trustworthy. A recommended solution for this common problem has been Bayesian model estimation. Bayesian methods rely on user specified information from historical data or researcher intuition to more accurately estimate the parameters.  This package provides a user friendly interface for estimating test reliability.  Here, reliability is modeled as a beta distributed random variable with shape parameters alpha=true score variance and beta=error variance.

Contact information
-------------------
Author: Joshua Ray Tanzer
<p>Email: jtanzer@lifespan.org

Functions
---------
* __bcov__: Given an N by P data matrix, this function will estimate a variance covariance matrix based on a prior inverse Wishart distribution.  Priors can be specified, if none are specified default prior variances are sample variances and default prior covariances are zero.  A Gibbs sampling approach is implemented.  Convergence should be assessed before interpreting results. Alpha and beta prior reliability parameters also need to be specified.  See below for a discussion of prior selection.  Alpha can be interpreted as prior true score variance and beta as prior error variance for the test.   For additional details on Bayesian methods and Monte Carlo Markov Chain procedures, see Hoff, P. D. (2009). A first course in Bayesian statistical methods. New York: Springer.
* __bcor__: Given an N by P data matrix, this function will estimate a correlation matrix based on a prior inverse Wishart distribution.  Priors can be specified, if none are specified default prior variances are sample variances and default prior covariances are zero.  A Gibbs sampling approach is implemented.  Convergence should be assessed before interpreting results.  Alpha and beta prior reliability parameters also need to be specified.  See below for a discussion of prior selection.  Alpha can be interpreted as prior true score variance and beta as prior error variance for the test.   For additional details on Bayesian methods and Monte Carlo Markov Chain procedures, see Hoff, P. D. (2009). A first course in Bayesian statistical methods. New York: Springer.
* __bomega__: This function is intended as an add on feature for measurement model estimation by blavaan.  If a blavaan model is specified with k items loading onto a single factor, bomega will extract the loadings and error variances to estimate coefficient omega from the model object.  Then, given prior true score variance (alpha) and prior error variance (beta), coefficient omega internal consistency reliability is estimated as a beta distributed random variable.  For more information on blavaan, see https://faculty.missouri.edu/~merklee/blavaan/.
* __bomega_general__: This function estimates coefficient omega internal consistency reliability from user specified item loadings and error variances, in addition to prior estimates of true score variance (alpha) and prior error variance (beta).
* __brxx_Cor__: When X any Y variables are specified with this function, reliability is estimated based on their correlation.  The bcov function is used (see above).  True score variance is estimated as the covariance between measures and error variance is estimated as SD_x*SD_2_y-Cov(x,y).  These values, plus alpha prior true score variance and beta prior error variance are used to estimate reliability as a beta distributed random variable.
* __brxx_Cor_general__: Based on user specified correlation and standard deviation estimates, as well as priors for true score variance (alpha) and error variance (beta), reliability is estimated as a beta distributed random variable.
* __brxx_general__: The most flexible function in this package, all that needs to be specified are estimates of sample and prior true score variance, as well as sample and prior error variance.  Median and credible limits is provided based on a beta distribution for reliability.
* __brxx_ICC__: This is an add on function for blme output.  If reliability is estimated from intraclass correlation coefficient, reliability is estimated given prior values for alpha (true score variance) and beta (error variance).  A blmer fitted model object from the blme package is provided, and the within subject effect and residual variance is extracted to estimate the ICC.  For more information on the blme package, see https://cran.r-project.org/web/packages/blme/blme.pdf.
* __brxx_ICC_general__: A general form for estimating reliability from intraclass correlation coefficient, all that needs to be specified are the within subject variance, residual variance, and prior shape parameters for reliability (alpha as true score variance and beta as error variance).
* __standardize__: This function standardizes all items in the dataset.  This is included in the package for ease of use.  It is strongly recommended to standardize items before combining scores into composites and estimating reliability, so as to avoid scaling effects on reliability estimates.  More discussion of this is provided below.
* __scree__: Calculates eigenvalues, percentage of variance, and their accumulative values.  Also provides the scree plot.  Eigenvalues are calculated based on the pairwise complete cases of the correlation matrix, to account for the possibility of missing data.
* __prep__: Creates a list object from a dataset for estimating a Bayesian exploratory factor analysis model.  All that needs to be specified is the data and the number of factors to be extracted.  Also able to be specified is the prior loading matrix for the means of the loadings (assumed to have variance IG(0.0001,0.0001)).  If no prior is given, the maximum likelihood solution will be assumed.  Stan code for the exploratory factor analysis model used in conjunction with this function is provided in the examples.
* __unpack__: Takes raw Stan output from the model included in the examples and converts it into parceled matrices for each parameter of the factor analysis model (Lambda matrix of item loadings, x_i matrix of latent trait scores, tau vector of item means, and alpha vector of loading variance estimates).  All that needs to be provided is an S by theta matrix of raw Stan output (see example) and list object of model information as is created by the format function.
* __process__: Provided with an unpacked loading matrix, this function calculates communality, uniqueness, and reliability estimates.  If the number of factors is greater than one, matrix rotation is performed and returns rotated loadings, interfactor correlations, loadings on the G factor, and estimates of hierarchical omega.  If the number of factors is one, total omega is calculated instead.  The S by P*Q loading matrix samples of the posterior needs to be provided, as well as the model information in list form from the format function.  Lastly, the matrix rotation can be specified based on the options in the GPArotation package; if roation is left blank, an oblimin rotation is assumed.  Outputs come as a list of S by Theta matrices for loadings, communalities, uniquenesses, and reliability estimates; if more than one factor is extracted, the list also includes G factor loadings and interfactor correlations.
* __summarize__: This function converts an S by Theta matrix of raw Stan output into summary matrices.  The raw data matrix of posterior samples must be specified, as well as the number of rows and columns for the target matrix.  This can come from the unpack statment (e.g. to summarize the x_i matrix of N by Q dimensionality of latent traits) or the process statement (e.g. to summarize the Lambda matrix of P by Q dimensionality of rotated item loadings).  An additional input is the CI, which specifies the width of credible intervals as a percentage.  If no limit is specified, 95% is assumed.  The output includes the median, standard deviation, upper limit and lower limit.  Limits on the CIs are highest posterior density based.  Lastly, a simple summary table is created, with the median (SD) and a * to indicate which values at the specified alpha level CI exclude a value of 0.

Prior selection
----------------
In this package, reliability is modeled as a beta distributed random variable.  This distribution was selected because of the clear concordance between how reliability is defined and the probable values for a beta distribution.  Scale reliability is defined as true score variance/(true score variance+error variance), and the range of possible values are from 0 to 1.  A beta distribution is defined by two parameters, alpha and beta, the range is also from 0 to 1, and the expected value is alpha/(alpha+beta).  It follows that modeling scale reliability as a beta distribution, alpha is interpreted as true score variance and beta as error variance.  This approach is more robust to small samples and missing at random data.  Further, the credible interval will be more appropriately scaled to the range of probable values; compared to a confidence interval, which assumes a normal distribution of errors.
<p>While this provides a methodical framework to more carefully model reliability, it is necessary to specify prior shape parameters, which may not be intuitive.  Historical data can be used to inform priors, however not all methods for estimating reliability clearly dilineate between true score variance and error variance (ie, reliability estimated from correlation).  As a solution, a more intuitive approach is proposed from beta quantiles.  If it is believed with 90% confidence that reliability will likely be between 0.40 and 0.90, then the shape parameters for the beta distribution containing quantile limits at the specified values can be solved for.  As a more intuitive solution, it is proposed that the prior could be based on the shape parameters for a beta distribution given specified quantile limits.  The [beta.parms.from.quantiles](http://www.medicine.mcgill.ca/epidemiology/joseph/pbelisle/R/BetaParmsFromQuantiles.R) function by Lawrence Joseph and Patrick Belisle is provided in the example to solve for these shape parameter values, to simplify prior selection estimation.

Item Scaling
------------
Because this method uses a beta distribution to infer the probable values of reliability, the larger the variance, true and error, the more precise the estimate.  It follows that two tests scaled at different ranges (e.g. 1 to 7 versus 1 to 100) would have vastly different reliability estimates as a function of their scaling minimum and maximum.  To address this, it is strongly recommended that each item be standardized before being summed to provide a composite score.  This will result in precision of estimate that is a function of test length rather than item scaling, which is already an expected and documented relationship.  The standardize function is provided to easily standardize item scores.  Additionally, the the user can specify the number of items on the test in the event that scores are already summed into composites.  If this is the case (e.g. if test retest reliability is desired but the individual item responses are missing), the composite scores can be standardized and the number of items on the test can be specified.  This will properly rescale the precision of reliability estimates to the test length.

Examples
--------
## Prior from quantiles
The [beta.parms.from.quantiles](http://www.medicine.mcgill.ca/epidemiology/joseph/pbelisle/R/BetaParmsFromQuantiles.R) function can be loaded as follows:

```{r}
beta.parms.from.quantiles <- function(q, p=c(0.025,0.975),
  precision=0.001, derivative.epsilon=1e-3, start.with.normal.approx=T, start=c(1, 1), plot=F)
{
  # Version 1.3 (February 2017)
  #
  # Function developed by 
  # Lawrence Joseph and pbelisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@rimuhc.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html
  #
  # Please refer to our webpage for details on each argument.
  
  f <- function(x, theta){dbeta(x, shape1=theta[1], shape2=theta[2])}
  F.inv <- function(x, theta){qbeta(x, shape1=theta[1], shape2=theta[2])}
  f.cum <- function(x, theta){pbeta(x, shape1=theta[1], shape2=theta[2])}
  f.mode <- function(theta){a <- theta[1]; b <- theta[2]; mode <- ifelse(a>1, (a-1)/(a+b-2), NA); mode}
  theta.from.moments <- function(m, v){a <- m*m*(1-m)/v-m; b <- a*(1/m-1); c(a, b)}
  plot.xlim <- c(0, 1)
  
  dens.label <- 'dbeta'
  parms.names <- c('a', 'b')
  
  if (length(p) != 2) stop("Vector of probabilities p must be of length 2.")
  if (length(q) != 2) stop("Vector of quantiles q must be of length 2.")
  p <- sort(p); q <- sort(q)
  
  #_____________________________________________________________________________________________________
        
  print.area.text <- function(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim, M=30, M0=50)
  {
    par.usr <- par('usr')
    par.din <- par('din')

    p.string <- as.character(round(c(0,1) + c(1,-1)*p.check, digits=4))
    str.width <- strwidth(p.string, cex=cex)
    str.height <- strheight("0", cex=cex)

    J <- matrix(1, nrow=M0, ncol=1)
    
    x.units.1in <- diff(par.usr[c(1,2)])/par.din[1]
    y.units.1in <- diff(par.usr[c(3,4)])/par.din[2]
    aspect.ratio <- y.units.1in/x.units.1in

    # --- left area  -----------------------------------------------------------

    scatter.xlim <- c(max(plot.xlim[1], par.usr[1]), q[1])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))

    grid.df <- f(x, theta)

    # Estimate mass center
    tmp.p <- seq(from=0, to=p[1], length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-1]*diff(tmp.x))/diff(range(tmp.x)))

    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}

    # Eliminate points to the right of the mode, if any
    w <- which(x>mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}

    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[1]+str.width[1]) <= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]

    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y, fromLast=T))
    if (length(w)){x <- x[-w]; y <- y[-w]}

    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[1], adj=c(1,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[1], mean(par.usr[c(3,4)]), labels=p.string[1], col='gray', cex=cex, srt=90, adj=c(1,0))
    }

    # --- right area  ----------------------------------------------------------

    scatter.xlim <- c(q[2], plot.xlim[2])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    grid.df <- f(x, theta)

    # Estimate mass center
    tmp.p <- seq(from=p[2], to=f.cum(plot.xlim[2], theta), length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-length(h)]*diff(tmp.x))/diff(range(tmp.x)))

    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}

    # Eliminate points to the left of the mode, if any
    w <- which(x<mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}

    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[2]-str.width[2]) >= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]

    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y))
    if (length(w)){x <- x[-w]; y <- y[-w]}

    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[2], adj=c(0,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[2], mean(par.usr[c(3,4)]), labels=p.string[2], col='gray', cex=cex, srt=-90, adj=c(1,0))
    }
  }
  
  # ......................................................................................................................................
  
  Newton.Raphson <- function(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start)
  {
    Hessian <- matrix(NA, 2, 2)

    if (start.with.normal.approx)
    {
      # Probably not a very good universal choice, but proved good in most cases in practice
      m <-  diff(q)/diff(p)*(0.5-p[1]) + q[1]
      v <- (diff(q)/diff(qnorm(p)))^2
      theta <- theta.from.moments(m, v)
    }
    else theta <- start
  
  
    change <- precision + 1
    niter <- 0
    # Newton-Raphson multivariate algorithm
    while (max(abs(change)) > precision)
    {
      Hessian[,1] <- (f.cum(q, theta) - f.cum(q, theta - c(derivative.epsilon, 0))) / derivative.epsilon
      Hessian[,2] <- (f.cum(q, theta) - f.cum(q, theta - c(0, derivative.epsilon))) / derivative.epsilon
  
      f <- f.cum(q, theta) - p
      change <- solve(Hessian) %*% f
      last.theta <- theta
      theta <- last.theta - change
  
      # If we step out of limits, reduce change
  
      if (any(theta<0))
      {
        w <- which(theta<0)
        k <- min(last.theta[w]/change[w])
        theta <- last.theta - k/2*change
      }
      
      niter <- niter + 1
    }
    
    list(theta=as.vector(theta), niter=niter, last.change=as.vector(change))
  }
  
  # ...............................................................................................................
  
  plot.density <- function(p, q, f, f.cum, F.inv, mode, theta, plot.xlim, dens.label, parms.names, cex)
  {
    if (length(plot.xlim) == 0)
    {
      plot.xlim <- F.inv(c(0, 1), theta)
      
      if (is.infinite(plot.xlim[1]))
      {
        tmp <- min(c(0.001, p[1]/10))
        plot.xlim[1] <- F.inv(tmp, theta)
      }  
      
      if (is.infinite(plot.xlim[2]))
      {
        tmp <- max(c(0.999, 1 - (1-p[2])/10))
        plot.xlim[2] <- F.inv(tmp, theta)
      }
    }
    plot.xlim <- sort(plot.xlim)
  
  
    x <- seq(from=min(plot.xlim), to=max(plot.xlim), length=1000)
    h <- f(x, theta)
    x0 <- x; f0 <- h
    ylab <- paste(c(dens.label, '(x, ', parms.names[1], ' = ', round(theta[1], digits=5), ', ', parms.names[2], ' = ', round(theta[2], digits=5), ')'), collapse='')
    plot(x, h, type='l', ylab=ylab)

    # fill in area on the left side of the distribution
    x <- seq(from=plot.xlim[1], to=q[1], length=1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from=max(plot.xlim), to=q[2], length=1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # draw distrn again
    points(x0, f0, type='l')
    h <- f(q, theta)
    points(rep(q[1], 2), c(0, h[1]), type='l', col='orange')
    points(rep(q[2], 2), c(0, h[2]), type='l', col='orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)  

    xaxp <- par("xaxp")
    x.ticks <- seq(from=xaxp[1], to=xaxp[2], length=xaxp[3]+1)
    q2print <- as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    mtext(q2print, side=1, col='orange', at=q2print, cex=0.6, line=2.1)
    points(q, rep(par('usr')[3]+0.15*par('cxy')[2], 2), pch=17, col='orange')
  }
  
  #________________________________________________________________________________________________________________


  parms <- Newton.Raphson(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start=start)
  p.check <- f.cum(q, parms$theta)
  
  if (plot) plot.density(p, q, f, f.cum, F.inv, f.mode(parms$theta), parms$theta, plot.xlim, dens.label, parms.names, 0.8)

  list(a=parms$theta[1], b=parms$theta[2], last.change=parms$last.change, niter=parms$niter, q=q, p=p, p.check=p.check)
}
```
So if you wanted 90% quantile limits (ie limits at 5% and 95%) that land at 0.4 and 0.9, then the following code can be used.
```{r}
CL=c(0.05,0.95)
Values=c(0.4,0.9)
beta.parms.from.quantiles(Values,CL)
```
The resultant values are alpha=3.51 and beta=1.75.

## bomega function with blavaan file
First, a measurement model needs to be fit using blavaan.
```{r}
mod='f=~X1+X2+X3+X4+X5'

fit=bsem(mod,data=your_data)
```
Then, read in the fitted blavaan model to the bomega statement. Here, K=5 because there are five items.
```{r}
bomega(mod=fit,k=5,alpha=3.51,beta=1.75,CI=0.95)
```

## brxx_ICC function from blme file
First, the ICC needs to be estimated using blme.
```{r}
fit=blmer(Y~(1|ID),data=your_data)
```
Then, read in the fitted blme model to the brxx_ICC statement.
```{r}
brxx_ICC(mod=fit,alpha=3.51,beta=1.75,CI=0.95)
```

## Bayesian factor analysis using Stan
To estimate a Bayesian factor analysis model, stan software can be used based on the following model code:
```{r}
fa="
data {
int<lower=0> N; // Number of observations
int<lower=0> P; // Number of items
int<lower=0> Q; // Number of latent dimensions
matrix[N, P] D_O; // Observed data matrix
matrix[P, Q] Load_Prior; // Prior loading matrix
matrix[P, P] I; //Identity matrix
matrix[N,P] R; //Missing data matrix
}

parameters {
matrix[N, Q] X; // Latent trait matrix
matrix[P, Q] Lambda; // Loading matrix
vector[P] tau; // Grand mean vector
vector<lower=0>[Q] alpha; // Loading variance vector
}

model {
matrix [N,P] X_Loading;
matrix [N,P] D_I;
for (i in 1:N) for (p in 1:P) 
X_Loading[i,p]=X[i,]*Lambda[p,]';
for (i in 1:N) for (p in 1:P) 
D_I[i,p]=D_O[i,p]*(1-R[i,p])+X_Loading[i,p]'*(R[i,p])+tau[p];

for (i in 1:N) for (q in 1:Q) X[i,q] ~ normal(0,1);
for (i in 1:N) tau ~ normal(0,1);
for(q in 1:Q) alpha[q] ~ inv_gamma(1e-4,1e-4);
for(q in 1:Q) Lambda[,q] ~ multi_normal(Load_Prior[,q], sqrt(alpha[q])*I);
for(i in 1:N) for (p in 1:P) D_I[i,p]~ normal(X_Loading[i,p]+tau[p], 1);
}

"

model=stan_model(model_code = fa)
```
This step may take some time, as Stan needs to rewrite the model code in C++.  Note that this model accounts for missing data by imputing the missing values with the expected values for individuals at each subsequent iteration of posterior sampling.  Raw data files can be prepared this model by the standardize and format functions.
```{R}
your_data_miss_s=standardize(your_data_miss)
formatted_data=prep(your_data_miss_s,nfactors=3)
```
Sampling the posterior distributions using this model can be saved in an S by Theta matrix using the following code.  This may also take some time, as the model reiterates parameter estimation.
```{R}
out=sampling(model, data=formatted_data, iter=5000, seed=999)
res=as.matrix(out)
```
Model convergence diagnostics can be summarized with
```{R}
plot(effectiveSize(Loading_Matrix),main="Effective Sample Size",ylab="N")
lines(rep(nrow(res),ncol(res)),col="red",lty=2)

plot(geweke.diag(res)$z,main="Geweke Diagnostic",ylab="Z")
lines(rep(1.96,ncol(res)),col="red",lty=2)
lines(rep(-1.96,ncol(res)),col="red",lty=2)
```
and individual trace plots can be examined using:
```{R}
plot(res[,1])
lines(res[,1])
```
All parameters should be visually examined for convergence.  The individual matrices of the factor analysis model can be separated out with the unpack function:
```{R}
unpacked=unpack(Samples=res,Format=formatted_data)
```
and the loadings can be rotated and further processed using the process function.
```{R}
processed=process(Loading_Matrix=unpacked$Loading_Matrix,
                  Format=formatted_data,
                  Rotate="oblimin")
```
Lastly, final summary tables of all relevant information can be created using
```{R}
summarize(processed$Loadings,
          nrow=Formatted_data$P,
          ncol=Formatted_data$Q)$Table
summarize(processed$Communality,
          nrow=Formatted_data$P,
          ncol=1)$Table
summarize(processed$Uniqueness,
          nrow=Formatted_data$P,
          ncol=1)$Table
summarize(processed$G_Factor,
          nrow=Formatted_data$P,
          ncol=1)$Table
summarize(processed$Interfactor_Correlations,
          nrow=Formatted_data$Q,
          ncol=Formatted_data$Q)$Table
summarize(processed$Omega,
          nrow=1,
          ncol=1)$Table
summarize(unpacked$Tau_Matrix,
          nrow=Formatted_data$P,
          ncol=1)$Table
```

Acknowledgements
----------------
Many thanks to Lisa Harlow for all her support and guidance.  Thanks to David Rindskopf for his recommendation on prior selection.  Thanks to Lawrence Joseph and Patrick Belisle for the beta.parms.from.quantiles function.  Lastly, thanks to Jing Wu for her guidance on Bayesian methods.
<p>Many ideas grow better when transplanted into another mind than the one where they sprang up

License
-------
[MIT](https://choosealicense.com/licenses/mit/)

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/JoshuaRayTanzer/brxx.svg?branch=master)](https://travis-ci.com/JoshuaRayTanzer/brxx)
<!-- badges: end -->
