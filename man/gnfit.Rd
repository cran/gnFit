\name{gnfit}
\alias{gnfit}
\title{Goodness of Fit Test for Continuous Distribution Functions}
\description{
Computes the test statistic and p-value of the Cramer-von Mises and Anderson-Darling test for some continuous distribution functions proposed by Chen and Balakrishnan (1995). In addition to our classic distribution functions here, we calculate the Goodness of Fit (GoF) test to dataset which follows the extreme value distribution function, without remembering the formula of distribution/density functions.
}

\usage{gnfit(dat, dist, df = 3, pr = NULL, threshold = NULL)}
\arguments{
   \item{dat}{A numeric vector of data values.}
   \item{dist}{A named distribution function in R, such as "norm", "t", "laplace", "logis", "gev", "gum", "gpd".}
    \item{df}{Degrees of freedom (> 2) for Student-t's distribution. This value is set to 3 by default.}
    \item{pr}{An object returned by maximum likelihood estimation of \href{https://CRAN.R-project.org/package=ismev}{gev.fit}, \href{https://CRAN.R-project.org/package=ismev}{gum.fit} or \href{https://CRAN.R-project.org/package=ismev}{gpd.fit}. It is also numeric vector giving the maximum likelihood estimation for the nonstationary model with location, scale and shape parameters for "gev" (generalized extreme value distribution), location and scale parameters in gumbel distribution or scale and shape parameters in case of "gpd" (generalized pareto distribution), resp.}
    \item{threshold}{The threshold is a single number. It is allocated by "gpd" with shape and scale parameters.}
}
\details{To test \eqn{H_0: X_1,\ldots, X_n} is a random sample from a continuous distribution with cumulative distribution function \eqn{F(x;\theta)}, where the form of \eqn{F} is known but \eqn{\theta} is unknown. We first esitmate \eqn{\theta} by \eqn{\theta^*} (for eg. maximum likelihood estimation method). Next, we compute \eqn{v_i=F(x_i,\theta^*)}, where the \eqn{x_i}'s are in ascending order.

The Cramer-von Mises test statistic is
\deqn{W^2=\sum_{i=1}^{n}(u_i-\frac{(2i-1)}{2n})^2+\frac{1}{12n},}
and
\deqn{A^2=-n-\frac{1}{n}\sum_{i=1}^{n}((2i-1)\log(u_i)+(2n+1-2i)\log(1-u_i));}
where, \eqn{u_i=\Phi((y_i-\bar{y}_i)/s_y)}, \eqn{y_i=\Phi^{-1}(v_i)}, and \eqn{\Phi} is the standard normal CDF and \eqn{\Phi^{-1}} its inverse.

Modify \eqn{W^2} and \eqn{A^2} into \deqn{W^*=W^2(1+0.5/n),} and  \deqn{A^*=A^2(1+0.75/n+2.25/n^2).}

The p-value is computed from the modified statistic \eqn{W^*} and \eqn{A^*} according to Table 4.9 in Stephens (1986).}

\value{The output is an object of the class "htest" for the Cramer-von Mises and Anderson-Darling statistics corresponding to p-values.}

\references{
Stephens (1986, ISBN:0824774876)

Chen and Balakrishnan (1995) <http://asq.org/qic/display-item/index.html?item=11407>

Marsaglia (2004) <doi:10.18637/jss.v009.i02>
}

\seealso{
The package of \href{https://CRAN.R-project.org/package=normtest}{nortest} for performing the Anderson-Darling test for normality. \href{ https://CRAN.R-project.org/package=ADGofTest}{ADGofTest} implementation of the Anderson-Darling goodness of fit test based on  Marsaglia's (2004).
}
\examples{
library(rmutil)
r <- rlaplace(1000, m = 1, s = 2)
gnfit(r, "laplace")

library(ismev)
pr <- c(-0.5, 1, 0.2)
r <- gevq(pr, runif(1000, 0, 1))
model <- gev.fit(r)$mle
gnfit(r, "gev", pr = model)
\donttest{
library(ismev)
r <- gum.q(runif(1000, 0, 1), -0.5, 1)
n <- length(r)
time <- matrix(1:n, ncol=1)
model <- gum.fit(r, ydat=time, mul=1)$mle
mle<-dim(2)
mle[1] <- model[1] + model[2] * (n+1)
mle[2] <- model[3]
gnfit(r, "gum", pr = mle)
}
\dontshow{
r <- rt (1000, df=4)
gnfit(r, "t", df=4)
}
}
