JDEoptim <- function(lower, upper, fn, constr = NULL, meq = 0, eps = 1e-5,
                     NP = 10*length(lower), Fl = 0.1, Fu = 1,
                     tau_F = 0.1, tau_CR = 0.1, tau_pF = 0.1,
                     jitter_factor = 0.001,
                     tol = 1e-15, maxiter = 2000*length(lower), fnscale = 1,
                     compare_to = c("median", "max"),
                     add_to_init_pop = NULL,
                     trace = FALSE, triter = 1,
                     details = FALSE, ...)

#   Copyright 2013, 2014, 2016, 2023, 2025, Eduardo L. T. Conceicao
#   Available under the GPL (>= 2)

{
    handle.bounds <- function(x, u) {
        # Check feasibility of bounds and enforce parameters limits
        # by a deterministic variant of bounce-back resetting
        # (also known as midpoint target/base)
        # Price, KV, Storn, RM, and Lampinen, JA (2005)
        # Differential Evolution: A Practical Approach to Global Optimization.
        # Springer, p 206
        bad <- x > upper
        x[bad] <- 0.5*(upper[bad] + u[bad])
        bad <- x < lower
        x[bad] <- 0.5*(lower[bad] + u[bad])
        x
    }

    performReproduction <- function() { # Mutate/recombine
        ignore <- runif(d) > CRtrial
        if (all(ignore))                  # ensure that trial gets at least
            ignore[sample(d, 1)] <- FALSE # one mutant parameter

        # Source for trial is the base vector plus weighted differential
        trial <- if (runif(1) <= pFtrial)
            X.base + Ftrial*(X.r1 - X.r2)
        else X.base + 0.5*(Ftrial + 1)*(X.r1 + X.r2 - 2*X.base)

        # or trial parameter comes from target vector X.i itself.
        trial[ignore] <- X.i[ignore]
        trial
    }

    which.best <- if (!is.null(constr))
        function(x) {
            ind <- TAVpop <= mu
            if (all(ind))
                which.min(x)
            else if (any(ind))
                which(ind)[which.min(x[ind])]
            else which.min(TAVpop)
        }
    else which.min

    # Check input parameters
    compare_to <- match.arg(compare_to)
    stopifnot(length(upper) == length(lower),
              length(lower) > 0, is.numeric(lower), is.finite(lower),
              length(upper) > 0, is.numeric(upper), is.finite(upper),
              lower <= upper,
              is.function(fn))
    if (!is.null(constr))
        stopifnot(is.function(constr),
                  length(meq) == 1, meq == as.integer(meq), meq >= 0,
                  is.numeric(eps), is.finite(eps), eps > 0,
                  length(eps) == 1 || length(eps) == meq)
    stopifnot(length(NP) == 1, NP == as.integer(NP), NP >= 0,
              length(Fl) == 1, is.numeric(Fl),
              length(Fu) == 1, is.numeric(Fu),
              Fl <= Fu)
    stopifnot(length(tau_F) == 1, is.numeric(tau_F), 0 <= tau_F, tau_F <= 1,
              length(tau_CR) == 1, is.numeric(tau_CR), 0 <= tau_CR, tau_CR <= 1,
              length(tau_pF) == 1, is.numeric(tau_pF), 0 <= tau_pF, tau_pF <= 1)
    if (!is.null(jitter_factor))
        stopifnot(length(jitter_factor) == 1,
                  is.numeric(jitter_factor),
                  is.finite(jitter_factor))
    stopifnot(length(tol) == 1, is.numeric(tol), is.finite(tol),
              length(maxiter) == 1, maxiter == as.integer(maxiter),
              maxiter >= 0,
              length(fnscale) == 1, is.numeric(fnscale),
              is.finite(fnscale), fnscale > 0)
    if (!is.null(add_to_init_pop))
        stopifnot(NROW(add_to_init_pop) == length(lower),
                  is.numeric(add_to_init_pop),
                  is.finite(add_to_init_pop),
                  add_to_init_pop >= lower,
                  add_to_init_pop <= upper)
    stopifnot(length(trace) == 1, is.logical(trace), !is.na(trace),
              length(triter) == 1, triter == as.integer(triter), triter >= 1,
              length(details) == 1, is.logical(details), !is.na(details))

    child <- if (is.null(constr)) { # Evaluate/select
        expression({
            ftrial <- fn1(trial)
            if (ftrial <= fpop[i]) {
                pop[, i] <- trial
                fpop[i] <- ftrial
                F[, i] <- Ftrial
                CR[i] <- CRtrial
                pF[i] <- pFtrial
            }
        })
    } else if (meq > 0) { # equality constraints are present
                          # alongside the inequalities
        # Zhang, Haibo, and Rangaiah, G. P. (2012).
        # An efficient constraint handling method with integrated differential
        # evolution for numerical and engineering optimization.
        # Computers and Chemical Engineering 37, 74-88.
        expression({
            htrial <- constr1(trial)
            TAVtrial <- sum( pmax(htrial, 0) )
            if (TAVtrial > mu) {
                if (TAVtrial <= TAVpop[i]) { # trial and target are both
                    pop[, i] <- trial        # infeasible, the one with smaller
                    hpop[, i] <- htrial      # constraint violation is chosen
                    F[, i] <- Ftrial         # or trial vector when both are
                    CR[i] <- CRtrial         # solutions of equal quality
                    pF[i] <- pFtrial
                    TAVpop[i] <- TAVtrial
                }
            } else if (TAVpop[i] > mu) { # trial is feasible and target is not
                pop[, i] <- trial
                fpop[i] <- fn1(trial)
                hpop[, i] <- htrial
                F[, i] <- Ftrial
                CR[i] <- CRtrial
                pF[i] <- pFtrial
                TAVpop[i] <- TAVtrial
            } else {                     # between two feasible solutions, the
                ftrial <- fn1(trial)     # one with better objective function
                if (ftrial <= fpop[i]) { # value is chosen
                    pop[, i] <- trial    # or trial vector when both are
                    fpop[i] <- ftrial    # solutions of equal quality
                    hpop[, i] <- htrial
                    F[, i] <- Ftrial
                    CR[i] <- CRtrial
                    pF[i] <- pFtrial
                    TAVpop[i] <- TAVtrial
                    FF <- sum(TAVpop <= mu)/NP
                    mu <- mu*(1 - FF/NP)
                }
            }
        })
    } else {              # only inequality constraints are present
        expression({
            htrial <- constr1(trial)
            TAVtrial <- sum( pmax(htrial, 0) )
            if (TAVtrial > mu) {
                if (TAVtrial <= TAVpop[i]) { # trial and target both infeasible
                    pop[, i] <- trial
                    hpop[, i] <- htrial
                    F[, i] <- Ftrial
                    CR[i] <- CRtrial
                    pF[i] <- pFtrial
                    TAVpop[i] <- TAVtrial
                }
            } else if (TAVpop[i] > mu) { # trial is feasible and target is not
                pop[, i] <- trial
                fpop[i] <- fn1(trial)
                hpop[, i] <- htrial
                F[, i] <- Ftrial
                CR[i] <- CRtrial
                pF[i] <- pFtrial
                TAVpop[i] <- TAVtrial
                FF <- sum(TAVpop <= mu)/NP
                mu <- mu*(1 - FF/NP)
            } else {                     # two feasible solutions
                ftrial <- fn1(trial)
                if (ftrial <= fpop[i]) {
                    pop[, i] <- trial
                    fpop[i] <- ftrial
                    hpop[, i] <- htrial
                    F[, i] <- Ftrial
                    CR[i] <- CRtrial
                    pF[i] <- pFtrial
                    TAVpop[i] <- TAVtrial
                    FF <- sum(TAVpop <= mu)/NP
                    mu <- mu*(1 - FF/NP)
                }
            }
        })
    }

    fn1 <- function(par) fn(par, ...)

    if (!is.null(constr))
        constr1 <-
            if (meq > 0) {
                eqI <- 1:meq
                function(par) {
                    h <- constr(par, ...)
                    h[eqI] <- abs(h[eqI]) - eps
                    h
                }
            } else function(par) constr(par, ...)

    use.jitter <- !is.null(jitter_factor)

    # Zielinski, Karin, and Laur, Rainer (2008).
    # Stopping criteria for differential evolution in
    # constrained single-objective optimization.
    # In: U. K. Chakraborty (Ed.), Advances in Differential Evolution,
    # SCI 143, Springer-Verlag, pp 111-138
    conv <- expression(
      ( do.call(compare_to, list(fpop)) - fpop[x.best.ind] )/fnscale
    )

    # Initialization
    d <- length(lower)
    pop <- matrix(runif(NP*d, lower, upper), nrow = d)
    if (!is.null(add_to_init_pop)) {
        pop <- unname(cbind(pop, add_to_init_pop))
        NP <- ncol(pop)
    }
    stopifnot(NP >= 4)
    # Combine jitter with dither
    # Storn, Rainer (2008).
    # Differential evolution research - trends and open questions.
    # In: U. K. Chakraborty (Ed.), Advances in Differential Evolution,
    # SCI 143, Springer-Verlag, pp 11-12
    F <- if (use.jitter)
        (1 + jitter_factor*runif(d, -0.5, 0.5)) %o% runif(NP, Fl, Fu)
    else matrix(runif(NP, Fl, Fu), nrow = 1)
    CR <- runif(NP)
    pF <- runif(NP)
    fpop <- apply(pop, 2, fn1)
    stopifnot(is.vector(fpop), !anyNA(fpop), !is.nan(fpop), !is.logical(fpop))
    if (!is.null(constr)) {
        hpop <- apply(pop, 2, constr1)
        stopifnot(is.matrix(hpop) || is.vector(hpop),
                  !anyNA(hpop), !is.nan(hpop), !is.logical(hpop))
        if (is.vector(hpop)) dim(hpop) <- c(1, length(hpop))
        TAVpop <- apply( hpop, 2, function(x) sum(pmax(x, 0)) )
        mu <- median(TAVpop)
    }

    popIndex <- 1:NP
    x.best.ind <- which.best(fpop)
    converge <- eval(conv)
    rule <- if (!is.null(constr))
        expression(converge >= tol || any(hpop[, x.best.ind] > 0))
    else expression(converge >= tol)
    convergence <- 0
    iteration <- 0

    while (eval(rule)) { # Generation loop
        if (iteration >= maxiter) {
            warning("maximum number of iterations reached without convergence")
            convergence <- 1
            break
        }
        iteration <- iteration + 1

        for (i in popIndex) { # Start loop through population

            # Equalize the mean lifetime of all vectors
            # Price, KV, Storn, RM, and Lampinen, JA (2005)
            # Differential Evolution: A Practical Approach to
            # Global Optimization. Springer, p 284
            i <- ((iteration + i) %% NP) + 1

            # Self-adjusting parameter control scheme
            Ftrial <- if (runif(1) <= tau_F) {
                # Combine jitter with dither
                if (use.jitter)
                    runif(1, Fl, Fu) * (1 + jitter_factor*runif(d, -0.5, 0.5))
                else runif(1, Fl, Fu)
            } else F[, i]

            CRtrial <- if (runif(1) <= tau_CR) runif(1) else CR[i]
            pFtrial <- if (runif(1) <= tau_pF) runif(1) else pF[i]

            # DE/rand/1/either-or/bin
            X.i <- pop[, i]
            # Randomly pick 3 vectors all diferent from target vector
            r <- sample(popIndex[-i], 3)
            X.base <- pop[, r[1L]]
            X.r1 <- pop[, r[2L]]
            X.r2 <- pop[, r[3L]]

            trial <- handle.bounds(performReproduction(), X.base)

            eval(child)

            x.best.ind <- which.best(fpop)
        }

        converge <- eval(conv)
        if (trace && (iteration %% triter == 0))
            cat(iteration, ":", "<", converge, ">", "(", fpop[x.best.ind], ")",
                pop[, x.best.ind],
                if (!is.null(constr))
                    paste("{", which(hpop[, x.best.ind] > 0), "}"),
                fill = TRUE)
    }

    res <- list(par = pop[, x.best.ind],
                value = fpop[x.best.ind],
                iter = iteration,
                convergence = convergence)
    if (details) {
        res$poppar <- pop
        res$popcost <- fpop
    }
    res
}

## Not exported, and only used because CRAN checks must be faster
doExtras <- function() {
    interactive() || nzchar(Sys.getenv("R_DEoptimR_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}
