% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pt_def_extras.R
\name{balance_intercepts}
\alias{balance_intercepts}
\title{Balancing intercepts: update a definitions table to (approximately) control marginal means}
\usage{
balance_intercepts(def, new_means, omit, default = 0, round_digits = 3)
}
\arguments{
\item{def}{Definition data.table to be modified}

\item{new_means}{Named numeric vector. Names are variables whose formula to modify; numbers are desired means.}

\item{omit}{Character vector. Names of variables to ignore (see \code{default})}

\item{default}{Numeric. Desired marginal mean for all variables not mentioned in \code{new_means}. Default 0.}

\item{round_digits}{Integer, default 3. Number of digits to round coefficients (trust me, you want this - if you don't round you will get a gigantic number of digits and the table will be unreadible.)}
}
\value{
A data.table that is a data definition as in \code{simstudy}.
}
\description{
Balancing intercepts: update a definitions table to (approximately) control marginal means
}
\examples{
formulas = list(
  w1 = 'w1{t-1} + a{t-1}', #note no intercepts - these will be added automatically
  w2 = 'w1{t}   + w2{t-1} + a{t-1}',
  w3 = 'w1{t}   + w2{t}   + w3{t-1} + a{t-1}',
  a  = 'w1{t}   + cos(w2{t})   + I(w3{t}^2)  + u',
  y  = 'sin(w1{t})   + w2{t}*w3{t}   + a{t} + u'
)
dists =  c(
 w1 = 'normal',
 w2 = 'normal',
 w3 = 'normal',
 a  = 'binary',
 y  = 'normal'
)
variances = c(
  w1 = .1,
  w2 = .1,
  w3=  .1,
  a = 1, #this has no effect
  y = 0.1
)
def = pt_def(formulas, dists, links=default_link(dists),
       variances, periods = 3, u_coef = 1)

balance_intercepts(def,
                   new_means = c(a1_conditional =0.3, a2_conditional = 0.3, y0=1, y1=1.5, y2=-1),
                   omit=c('a1','a2'),
                   default = 0)
}
