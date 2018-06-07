\name{rskFac}
\alias{rskFac}
\title{Risk Factors}
\description{
The Value at Risk (VaR) of level \eqn{\alpha} (\eqn{\alpha}-quantile) of an event is a number attempting to summarize the risk of that event and define the worst expected loss of the event over a period of time. The Average VaR is another important measure of the risk at a given confidence level, which calculated by using the function of "rskFac".
}

\usage{rskFac(dat, alpha = 0.1, dist = "norm", df = NULL)}
\arguments{
   \item{dat}{A numeric vector of object data.}
   \item{alpha}{Confidence level \eqn{\alpha} (\eqn{0<\alpha<0.5}).}
   \item{dist}{A named of distribution function which should be fitted to data values. The distibution function is selected by the name of "laplace", "logis", "gum", "t" and "norm".}
   \item{df}{degrees of freedom from a specified distribution function.}
}
\details{Suppose \eqn{X} is random variable (rv) has distribution function (df) \eqn{F}. Given a confidence level \eqn{\alpha\in (0, 1),} Value at Risk (VaR) of the underlying \eqn{X} at the confidence level \eqn{\alpha} is the smallest number \eqn{x}  such that the probability that the underlying \eqn{X} exceeds \eqn{x} is at least \eqn{1-\alpha.}
	In other word, if \eqn{X} is a rv with symmetric distribution function \eqn{F} (e.g., the return value of a portfolio), then \eqn{VaR_{\alpha}} is the negative of the \eqn{\alpha} quantile, i.e.,
\deqn{VaR_{\alpha}(X)=Q(\alpha)=inf{x \in Real : Pr( X \le  x )\le \alpha}.
}
	where, \eqn{Q(.)=F^{-1}(.).}

Since, the \eqn{VaR_\alpha(X)} is the nagative of \eqn{\alpha} quantile in the left tail, \eqn{-VaR_{1-\alpha}(-X)} is positive value of VaR in right tail.

	The average \eqn{VaR_\alpha,} \eqn{(AVaR_\alpha)} for \eqn{0<\alpha\le 1} of \eqn{X} is defined as
\deqn{
	AVaR_\alpha(X)= \frac{1}{\alpha}\int_{0}^{\alpha}VaR(x) dx,
}
	The AVaR is known under the names of conditional VaR (CVaR), tail VaR (TVaR) and expected shortfall.

	Pflug and Romisch (2007, ISBN: 9812707409) shows the AVaR may be represented as the optimal value of the following optimization problem
	\deqn{AVaR_\alpha (X) = VaR_\alpha(X) - \frac{1}{\alpha} E((X - VaR_\alpha(X))^{-}).}
	where, \eqn{(y)^{-} = min (y,0)}.
	To approximate the integral, it is given by
	\deqn{AVaR_\alpha(X)=VaR_\alpha(X)+\frac{1}{t \alpha}\sum_{i=1}^{t}max{(VaR_\alpha(X) - X), 0},}
	where, \eqn{t} is number of observations. By considering the rv \eqn{-X}, the \eqn{-AVaR_{1-\alpha}} in right tail is obtainable.
}
\value{
The values of output are "VaR", "AVaR_n" and "AVaR_p" correspond to the VaR, Average VaR in left tail, Average VaR in right tail.
}
\references{
Pflug and Romisch (2007, ISBN: 9812707409)

}

\examples{
library(rmutil)
r <- rlaplace(1000, m = 1, s = 2)
rskFac(r, dist = "laplace", alpha = 0.1)

}
