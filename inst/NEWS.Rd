\name{NEWS}
\title{News for package \pkg{DEoptimR}}

\section{Changes in DEoptimR version 1.1-0 (2023-07-10)}{
  \itemize{
    \item New \code{NCDEoptim()} adds a bespoke implementation of the NCDE
          algorithm for multimodal optimization.
  }
}

\section{Changes in DEoptimR version 1.0-6 (2016-07-05)}{
  \itemize{
    \item In \code{JDEoptim()}, \code{FUN} renamed to \code{compare_to}.
  }
}

\section{Changes in DEoptimR version 1.0-5 (2016-07-01)}{
  \itemize{
    \item In \code{JDEoptim()}, \code{tau1}, \code{tau2} and \code{tau3}
          renamed to \code{tau_F}, \code{tau_CR} and \code{tau_pF},
          respectively.
  }
}

\section{Changes in DEoptimR version 1.0-0 (2014-01-27)}{
  \itemize{
    \item First release of DEoptimR: implementation of a bespoke variant
          of the jDE algorithm for single-objective optimization.
    \item Constraint handling based on biasing feasible over unfeasible
          solutions by a parameter-less variant of the
          \eqn{\varepsilon}{epsilon}-constrained method.
  }
}