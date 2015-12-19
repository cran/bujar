\name{bujar}
\alias{bujar}
\alias{print.bujar}
\alias{plot.bujar}
\alias{summary.bujar}
\title{ Buckley-James Regression }
\description{
Buckley-James regression for right-censoring survival data with high-dimensional covariates. Including L_2 boosting with componentwise linear least squares, componentwise P-splines, regression trees. Other Buckley-James methods including elastic net, MCP, SCAD, MARS and ACOSSO (ACOSSO not supported for the current version). 
}
\usage{
bujar(y, cens, x, valdata = NULL, degree = 1, learner = "linear.regression",
center=TRUE, mimpu = NULL, iter.bj = 20, max.cycle = 5, nu = 0.1, mstop = 50, 
twin = FALSE, mstop2= 100, tuning = TRUE, cv = FALSE, nfold = 5, method = "corrected", 
vimpint = TRUE,gamma = 3, lambda=NULL, whichlambda=NULL, lamb = 0, s = 0.5, nk = 4, 
wt.pow = 1, theta = NULL, rel.inf = FALSE, tol = .Machine$double.eps, n.cores= 2, 
rng=123, trace = FALSE)
}
\arguments{
  \item{y}{ survival time}
  \item{cens}{ censoring indicator, must be 0 or 1 with 0=alive, 1=dead}
  \item{x}{ covariate matrix}
  \item{valdata}{ test data, which must have the first column as survival time, second column as censoring indicator, and the remaining columns similar to same x.}
  \item{degree}{ mars/tree/linear regression degree of interaction; if 2, second-order interaction,
                                          if degree=1, additive model;}
  \item{learner}{ methods used for BJ regression.} 
  \item{center}{ center covariates}
  \item{mimpu}{ initial estimate. If TRUE, mean-imputation;
  FALSE, imputed with the marginal best variable linear regression; if NULL, 0.}
  \item{iter.bj}{ number of B-J iteration}
  \item{max.cycle}{ max cycle allowed }
  \item{nu}{ step-size boosting parameter}
  \item{mstop}{ boosting tuning parameters. It can be one number or have the length \code{iter.bj}+\code{max.cycle}. If \code{cv=TRUE}, then \code{mstop} is the maximum number of tuning parameter}
  \item{twin}{ logical, if TRUE, twin boosting}
  \item{mstop2}{ twin boosting tuning parameter}
  \item{tuning}{ logical value. if TRUE, the tuning parameter will be selected by cv or AIC/BIC methods. Ignored if \code{twin=TRUE} for which no tuning parameter selection is implemented }
  \item{cv}{ logical value. if TRUE, cross-validation for tuning parameter, only used if \code{tuning=TRUE}. If \code{tuning=FALSE} or \code{twin=TRUE}, then ignored}
  \item{nfold}{ number of fold of cv}
  \item{method}{ boosting tuning parameter selection method in AIC}
  \item{vimpint}{ logical value. If TRUE, compute variable importance and interaction measures for MARS if \code{learner="mars"} and \code{degree} > 1. }
  \item{gamma}{ MCP, or SCAD gamma tuning parameter}
  \item{lambda}{ MCP, or SCAD lambda tuning parameter}
  \item{whichlambda}{ which lambda used for MCP or SCAD lambda tuning parameter }
  \item{lamb}{ elastic net lambda tuning parameter, only used if \code{learner="enet"}}
  \item{s}{ the second enet tuning parameter, which is a fraction between (0, 1), only used if \code{learne="enet"}}
  \item{nk}{ number of basis function for \code{learner="mars"}}
  \item{wt.pow}{ not used but kept for historical reasons, only for \code{learner=ACOSSO}. This is a parameter (power of weight). It might be chosen by CV from c(0, 1.0, 1.5, 2.0, 2.5, 3.0). If wt.pow=0, then this is COSSO method} 
  \item{theta}{ For \code{learner="acosso"}, not used now. A numerical vector with 0 or 1. 0 means the variable not included and 1 means included. See Storlie et al. (2009).}
  \item{rel.inf}{ logical value. if TRUE, variable importance measure and interaction importance measure computed}
  \item{tol}{ convergency criteria }
\item{n.cores}{The number of CPU cores to use. The cross-validation loop
  will attempt to send different CV folds off to different cores. Used for \code{learner="tree"}} 
  \item{rng}{ a number to be used for random number generation in boosting trees }
  \item{trace}{ logical value. If TRUE, print out interim computing results}
}
\details{
 Buckley-James regression for right-censoring survival data with high-dimensional covariates. Including L_2 boosting with componentwise linear least squares, componentwise P-splines, regression trees. Other Buckley-James methods including elastic net, SCAD and MCP. \code{learner="enet"} and \code{learner="enet2"} use two different implementations of LASSO. Some of these methods are discussed in Wang and Wang (2010) and the references therein. Also see the references below.
}
\value{
  \item{x}{ original covariates}
  \item{y}{ survival time}
  \item{cens}{ censoring indicator}
  \item{ynew}{ imputed y}
  \item{yhat}{ estimated y from ynew}
  \item{pred.bj}{ estimated y from the testing sample}
  \item{res.fit}{ model fitted with the learner}
  \item{learner}{ original learner used}
  \item{degree}{ =1, additive model, degree=2, second-order interaction}
  \item{mse}{ MSE at each BJ iteration, only available in simulations, or when valdata provided}
  \item{mse.bj}{ MSE from training data at the BJ termination}
  \item{mse.bj.val}{MSE with valdata}
  \item{mse.all}{ a vector of MSE for uncensoring data at BJ iteration}
  \item{nz.bj.iter}{ number of selected covariates at each BJ iteration}
  \item{nz.bj}{ number of selected covariates at the claimed BJ termination}
  \item{xselect}{ a vector of dimension of covariates, either 1 (covariate selected) or 0 (not selected)}
  \item{coef.bj}{ estimated coefficients with linear model}
  \item{vim}{ a vector of length of number of column of x, variable importance, between 0 to 100}
  \item{interactions}{ measure of strength of interactions}
  \item{ybstdiff}{ largest absolute difference of estimated y. Useful to monitor convergency}
  \item{ybstcon}{ a vector with length of BJ iteration each is a convergency measure} 
  \item{cycleperiod}{ number of cycle of BJ iteration}
  \item{cycle.coef.diff}{ within cycle of BJ, the maximum difference of coefficients for BJ boosting}
  \item{nonconv}{ logical value. if TRUE, non-convergency}
  \item{fnorm2}{ value of L_2 norm, can be useful to access convergencey}
  \item{mselect}{ a vector of length of BJ iteration, each element is the tuning parameter mstop}
  \item{contype}{ 0 (converged), 1, not converged but cycle found, 2, not converged and max iteration reached.}
}
\references{ 

  Zhu Wang and C.Y. Wang (2010),
  Buckley-James Boosting for Survival Analysis with High-Dimensional
Biomarker Data.
  \emph{Statistical Applications in Genetics and Molecular Biology},
  Vol. 9 : Iss. 1, Article 24.
 
  Peter Buhlmann and Bin Yu (2003),
  Boosting with the L2 loss: regression and classification.
  \emph{Journal of the American Statistical Association}, \bold{98},
  324--339.

  Peter Buhlmann (2006), Boosting for high-dimensional linear models.
  \emph{The Annals of Statistics}, \bold{34}(2), 559--583.

  Peter Buhlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, \bold{22}(4), 477--505.
  
  J. Friedman (1991), Multivariate Adaptive Regression Splines (with
  discussion) .
  \emph{Annals of Statistics}, \bold{19}/1, 1--141.

  J.H. Friedman, T. Hastie and R. Tibshirani (2000), Additive Logistic Regression:
a Statistical View of Boosting. \emph{Annals of Statistics} \bold{28}(2):337-374.

  C. Storlie, H. Bondell, B. Reich and H. H. Zhang (2009),
  Surface Estimation, Variable Selection, and the Nonparametric Oracle
        Property.
  \emph{Statistica Sinica}, to appear.

  Sijian Wang, Bin Nan, Ji Zhu, and David G. Beer (2008),
  Doubly penalized Buckley-James Method for Survival Data with High-Dimensional
        Covariates.
  \emph{Biometrics},
  \bold{64}:132-140.

  H. Zou and T. Hastie (2005), Regularization and variable selection via the elastic net.
  \emph{Journal of the Royal Statistical Society}, Series B, \bold{67}, 301-320.

}
\author{ Zhu Wang }
\examples{
data("wpbc", package = "TH.data")
wpbc2 <- wpbc[, 1:12]
wpbc2$status <- as.numeric(wpbc2$status) - 1
fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x= wpbc2[, -(1:2)])
print(fit)
coef(fit)
pr <- predict(fit)
plot(fit)
fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x= wpbc2[, -(1:2)], tuning = TRUE)
\dontrun{
fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x=wpbc2[, -(1:2)], learner="pspline")
fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x=wpbc2[, -(1:2)], 
 learner="tree", degree=2)
### select tuning parameter for "enet"
tmp <- gcv.enet(y=log(wpbc2$time), cens=wpbc2$status, x=wpbc2[, -(1:2)])
fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x=wpbc2[, -(1:2)], learner="enet", 
lamb = tmp$lambda, s=tmp$s)

fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x=wpbc2[, -(1:2)], learner="mars", 
degree=2)
summary(fit)
}
}
