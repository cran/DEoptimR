\name{NCDEoptim}
\alias{NCDEoptim}
\title{
  Bound-Constrained and Nonlinear Constrained Multimodal Optimization
  via Differential Evolution
}
\description{
  A bespoke implementation of the \sQuote{NCDE}
  (neighborhood based crowding DE) algorithm by
  Qu \emph{et al.} (2012) \doi{10.1109/TEVC.2011.2161873},
  assisted with the dynamic archive mechanism of
  Epitropakis \emph{et al.} (2013) \doi{10.1109/CEC.2013.6557556}.
}
\usage{
NCDEoptim(lower, upper, fn,
          constr = NULL, meq = 0, eps = 1e-5,
          crit = 1e-5, niche_radius = NULL, archive_size = 100,
          reinit_if_solu_in_arch = TRUE,
          NP = 100, Fl = 0.1, Fu = 1, CRl = 0, CRu = 1.1,
          nbngbrsl = NP/20, nbngbrsu = NP/5,
          tau_F = 0.1, tau_CR = 0.1, tau_pF = 0.1,
          tau_nbngbrs = 0.1,
          jitter_factor = 0.001,
          maxiter = 2000,
          add_to_init_pop = NULL,
          trace = FALSE, triter = 1,
          ...)
}
\arguments{
  \item{lower, upper}{numeric vectors, the lower and upper bounds of the
    search space (\emph{box constraints}); must be finite
    (\code{\link{is.finite}}).}

  \item{fn}{a \code{\link{function}} to be \strong{minimized} that takes a
    numeric vector \eqn{X_i} as first argument and returns the value of the
    objective.}

  \item{constr}{a vector \code{\link{function}} specifying the
    \strong{left-hand side} of equality constraints defined to equal zero
    (\eqn{h_j(X_i) = 0,\; j = 1,\ldots,\mathrm{meq}}),
    followed by inequality constraints defined as lesser than zero
    (\eqn{g_j(X_i) \le 0,\; j = \mathrm{meq}+1,\ldots}).  This function takes
    \eqn{X_i} as its first argument and returns a numeric vector with the same
    length of the total number of constraints.  It defaults to \code{NULL},
    which means that \strong{bound-constrained} minimization is used.}

  \item{meq}{an integer, the first \code{meq} constraints are
    \emph{equality} constraints whereas the remaining ones are
    \emph{inequality} constraints.  Defaults to \code{0}
    (inequality constraints only).}

  \item{eps}{the maximal admissible constraint violation for
    equality constraints.  A numeric vector of small positive tolerance values
    with length \code{meq} used in the transformation of equalities into
    inequalities of the form \eqn{|h_j(X_i)| - \epsilon \le 0}.  A scalar value
    is expanded to apply to all equality constraints.  Default is \code{1e-5}.}

  \item{crit}{a numeric, the acceptance threshold on the archive strategy.  If
    \code{\link{isTRUE}(\link{all.equal}(fn(X_best_so_far_in_archive), fn(X_i), tolerance = crit))},
    a solution \eqn{X_i} is checked for possible insertion into the
    dynamic archive.  Defaults to \code{1e-5}.}

  \item{niche_radius}{a numeric, the absolute tolerance used to decide whether
    the solution \eqn{X_i} is \emph{identical} to an already existing
    local or global solution \emph{in the archive}.  It defaults to \code{NULL},
    meaning that the niche radius is adaptively chosen during the search.
    Results are \strong{much better} if one is able to provide
    a reasonable value.}

  \item{archive_size}{an integer, the maximum number of solutions that
    can be kept in the archive; entries above this limit are discarded.
    Default is \code{100}.}

  \item{reinit_if_solu_in_arch}{a logical, if \code{TRUE}, any solution
    \eqn{X_i} already in the archive \strong{reinitializes} its nearest neighbor
    \emph{in the population} within the range
    \eqn{[\mathrm{lower}, \mathrm{upper}]}.  Default is \code{TRUE}.}

  \item{NP}{an integer, the population size.  Defaults to \code{100}.}

  \item{Fl}{a numeric, the minimum value that the
    \emph{scaling factor} \code{F} could take.  It defaults to \code{0.1}.}

  \item{Fu}{a numeric, the maximum value that the
    \emph{scaling factor} \code{F} could take.  It defaults to \code{1}.}

  \item{CRl}{a numeric, the minimum value to be used for the
    \emph{crossover constant} \code{CR}.  It defaults to \code{0}.}

  \item{CRu}{a numeric, the maximum value to be used for the
    \emph{crossover constant} \code{CR}.  It defaults to \code{1.1}.}

  \item{nbngbrsl}{an integer, the lower limit for the
    \emph{neighborhood size} \code{nbngbrs}.
    It defaults to \code{1/20} of the population size.}

  \item{nbngbrsu}{an integer, the upper limit for the
    \emph{neighborhood size} \code{nbngbrs}.
     It defaults to \code{1/5} of the population size.}

  \item{tau_F}{a numeric, the probability that the
    \emph{scaling factor} \code{F} is updated.  Defaults to \code{0.1}.}

  \item{tau_CR}{a numeric, the probability that the
    \emph{crossover constant} \code{CR} is updated.  Defaults to \code{0.1}.}

  \item{tau_pF}{a numeric, the probability that the
    \emph{mutation probability} \eqn{p_F}{pF} in the mutation strategy
    DE/rand/1/either-or is updated.  Defaults to \code{0.1}.}

  \item{tau_nbngbrs}{a numeric, the probability that the
    \emph{neighborhood size} \code{nbngbrs} is updated.  Defaults to \code{0.1}.}

  \item{jitter_factor}{a numeric, the tuning constant for \emph{jitter}.
    If \code{NULL} only \emph{dither} is used.  Default is \code{0.001}.}

  \item{maxiter}{an integer, the maximum number of iterations allowed which is
    the \strong{stopping condition}.  Default is \code{2000}.}

  \item{add_to_init_pop}{numeric vector of length \code{length(lower)} or
    column-wise \code{\link{matrix}} with \code{length(lower)} rows specifying
    initial candidate solutions which are appended to the randomly generated
    initial population.  Default is \code{NULL}.}

  \item{trace}{a logical, determines whether or not to monitor the
    iteration process.  Default is \code{FALSE}.}

  \item{triter}{an integer, trace output is printed at every
    \code{triter} iterations.  Default is \code{1}.}

  \item{\dots}{additional arguments passed to \code{fn} and \code{constr}.}
}
\details{
  This implementation differs mainly from the original \sQuote{NCDE} algorithm
  of Qu \emph{et al.} (2012) by employing the archiving procedure proposed
  in Epitropakis \emph{et al.} (2013) and the adaptive \sQuote{jDE} strategy
  instead of canonical Diferential Evolution.  The key reason for archiving
  good solutions during the search process is to prevent them from being lost
  during evolution.  Constraints are tackled through the
  \eqn{\varepsilon}{epsilon}-constrained method as proposed
  in Poole and Allen (2019).  The \sQuote{jDE} and
  \eqn{\varepsilon}{epsilon}-constrained mechanisms are applied in the same way
  as in \code{\link{JDEoptim}}, but with \emph{synchronous} mode of
  population update.  In contrast, the reinitialization in the current
  population triggered by already found solutions is done \emph{asynchronously}.

  Each line of trace output follows the format of:

  \code{iteration : < value of niche radius > population>> ( value of best solution ) best solution { index of violated constraints } archive>> [ number of solutions found ] ( value of best solution ) best solution}
}
\value{
  A list with the following components:
  \item{solution_arch}{a \code{\link{matrix}} whose columns are the
    local and global minima stored in the \strong{archive} of feasible solutions
    in ascending order of the objective function values.}

  \item{objective_arch}{the values of \eqn{\mathrm{fn}(X_i)} for the
    corresponding columns of \code{solution_arch}.}

  \item{solution_pop}{a \code{\link{matrix}} whose columns are the
    local and global minima stored in the \strong{final population}
    in ascending order of the objective function values;
    feasible solutions come first followed by the infeasible ones.}

  \item{objective_pop}{the values of \eqn{\mathrm{fn}(X_i)} for the
    corresponding columns of \code{solution_pop}.}

  \item{iter}{the number of iterations used.}
  and if there are general constraints present:
  \item{constr_value_arch}{a \code{\link{matrix}} whose columns contain
    the values of the constraints for \code{solution_arch}.}

  \item{constr_value_pop}{a \code{\link{matrix}} whose columns contain
    the values of the constraints for \code{solution_pop}.}
}
\references{
  Epitropakis, M. G., Li, X. and Burke, E. K. (2013)
  A dynamic archive niching differential evolution algorithm
  for multimodal optimization;
  in \emph{2013 IEEE Congress on Evolutionary Computation (CEC)}.
  IEEE, pp. 79--86.
  \doi{10.1109/CEC.2013.6557556}.

  Poole, D. J. and Allen, C. B. (2019)
  Constrained niching using differential evolution.
  \emph{Swarm and Evolutionary Computation} \bold{44}, 74--100.
  \doi{10.1016/j.swevo.2018.11.004}.

  Qu, B. Y., Suganthan, P. N. and Liang, J. J. (2012)
  Differential evolution with neighborhood mutation for multimodal optimization.
  \emph{IEEE Transactions on Evolutionary Computation} \bold{16}, 601--614.
  \doi{10.1109/TEVC.2011.2161873}.
}
\author{
  Eduardo L. T. Conceicao \email{mail@eduardoconceicao.org}
}
\note{
  \bold{This function is in an experimental stage.}
}
\examples{
\donttest{
# NOTE: Examples were excluded from testing
#       to reduce package check time.

# Use a preset seed so test values are reproducible.
set.seed(1234)

# Warning: the examples are run using a very small number of
# iterations to decrease execution time.

# Bound-constrained optimization

#   Vincent function
#
#   f(x) = -mean(sin(10*log(x)))
#
#   0.25 <= xi <= 10, i = {1, 2, ..., n}
#   The function has 6^n global minima without local minima.

NCDEoptim(c(0.25, 0.25), c(10, 10),
          function(x) -mean(sin(10*log(x))),
          niche_radius = 0.2,
          maxiter = 200, trace = TRUE, triter = 20)

# Nonlinear constrained optimization

#   Function F10 of Poole and Allen (2019)
#
#   f(x) = -sin(5*pi*x)^6 + 1
#   subject to:
#   g(x) = -cos(10*pi*x) <= 0
#
#   0 <= x <= 1
#   The 10 global optima are
#   (x1*, ..., x10*; f*) = ((2*(1:10) - 1)/20; 0.875).

NCDEoptim(0, 1,
          function(x) -sin(5*pi*x)^6 + 1,
          function(x) -cos(10*pi*x),
          niche_radius = 0.05,
          maxiter = 200, trace = TRUE, triter = 20)
}
}
\concept{multimodal optimization}
