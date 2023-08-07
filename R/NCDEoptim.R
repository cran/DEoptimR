NCDEoptim <- function(
    lower, upper, fn, constr = NULL, meq = 0, eps = 1e-5,
    crit = 1e-5, niche_radius = NULL, archive_size = 100,
    reinit_if_solu_in_arch = TRUE,
    NP = 100, Fl = 0.1, Fu = 1, CRl = 0, CRu = 1.1,
    nbngbrsl = NP/20, nbngbrsu = NP/5,
    tau_F = 0.1, tau_CR = 0.1, tau_pF = 0.1,
    tau_nbngbrs = 0.1,
    jitter_factor = 0.001,
    maxiter = 2000,
    add_to_init_pop = NULL, trace = FALSE, triter = 1,
    ...) {

#   Copyright 2023, Eduardo L. T. Conceicao
#   Available under the GPL (>= 2)

    handle_bounds <- function(x, u) {
        # Check feasibility of bounds and enforce parameters limits
        # by a deterministic variant of bounce-back resetting
        # Price, KV, Storn, RM, and Lampinen, JA (2005)
        # Differential Evolution: A Practical Approach to Global Optimization.
        # Springer, p 206
        bad <- x > upper
        x[bad] <- 0.5*(upper[bad] + u[bad])
        bad <- x < lower
        x[bad] <- 0.5*(lower[bad] + u[bad])
        x
    }

    perform_reproduction <- function() { # Mutate/recombine
        ignore <- runif(d) > CRtrial
        if (all(ignore))                  # ensure that trial gets at least
            ignore[sample(d, 1)] <- FALSE # one mutant parameter
        # Source for trial is the base vector plus weighted differential
        trial <- if (runif(1) <= pFtrial)
            X_base + Ftrial*(X_r1 - X_r2)
        else X_base + 0.5*(Ftrial + 1)*(X_r1 + X_r2 - 2*X_base)
        # or trial parameter comes from target vector X_i itself.
        trial[ignore] <- X_i[ignore]
        trial
    }

    which_best <- if (!is.null(constr))
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
    stopifnot(length(upper) == length(lower),
              is.numeric(lower), is.finite(lower),
              is.numeric(upper), is.finite(upper),
              lower <= upper,
              is.function(fn))
    if (!is.null(constr))
        stopifnot(is.function(constr),
                  length(meq) == 1, meq == as.integer(meq), meq >= 0,
                  is.numeric(eps), is.finite(eps), eps > 0,
                  length(eps) == 1 || length(eps) == meq)
    stopifnot(length(crit) == 1, is.numeric(crit),
              is.finite(crit), crit > 0)
    if (!is.null(niche_radius))
        stopifnot(length(niche_radius) == 1, is.numeric(niche_radius),
                  is.finite(niche_radius), niche_radius > 0)
    stopifnot(length(archive_size) == 1,
              archive_size == as.integer(archive_size),
              archive_size >= 0,
              length(reinit_if_solu_in_arch) == 1,
              is.logical(reinit_if_solu_in_arch))
    stopifnot(length(NP) == 1, NP == as.integer(NP), NP >= 0,
              length(Fl) == 1, is.numeric(Fl),
              length(Fu) == 1, is.numeric(Fu),
              Fl <= Fu,
              length(CRl) == 1, is.numeric(CRl),
              length(CRu) == 1, is.numeric(CRu),
              CRl <= CRu)
    stopifnot(length(tau_F) == 1, is.numeric(tau_F), tau_F >= 0, tau_F <= 1,
              length(tau_CR) == 1, is.numeric(tau_CR), tau_CR >= 0, tau_CR <= 1,
              length(tau_pF) == 1, is.numeric(tau_pF), tau_pF >= 0, tau_pF <= 1,
              length(tau_nbngbrs) == 1, is.numeric(tau_nbngbrs),
              tau_nbngbrs >= 0, tau_nbngbrs <= 1)
    if (!is.null(jitter_factor))
        stopifnot(length(jitter_factor) == 1,
                  is.numeric(jitter_factor),
                  is.finite(jitter_factor))
    stopifnot(length(maxiter) == 1,
              maxiter == as.integer(maxiter),
              maxiter >= 0)
    if (!is.null(add_to_init_pop))
        stopifnot(NROW(add_to_init_pop) == length(lower),
                  is.numeric(add_to_init_pop),
                  is.finite(add_to_init_pop),
                  add_to_init_pop >= lower,
                  add_to_init_pop <= upper)
    stopifnot(length(trace) == 1, is.logical(trace),
              length(triter) == 1, triter == as.integer(triter), triter >= 1)

    check_archive <- if (reinit_if_solu_in_arch) {
        expression({
            if (ftrial < best_fpop ||
                isTRUE(all.equal(best_fpop, ftrial, tolerance = crit))) {
                if (is.null(S))
                    S <- as.matrix(c(ftrial, trial))
                if (ftrial < best_fpop)
                    best_fpop <- ftrial
                found_ind <- sqrt(colSums(
                    (trial - S[-1, , drop = FALSE])^2
                )) <= R
                if (any(found_ind)) {
                    # Re-initialize nearest neighbor of the trial vector
                    pop[, k] <- runif(d, lower, upper)
                    fpop[k] <- fn1(pop[, k])
                    F[, k] <- if (use_jitter)
                        runif(1, Fl, Fu) * (1 + jitter_factor*runif(d, -0.5, 0.5))
                    else runif(1, Fl, Fu)
                    CR[k] <- runif(1, CRl, CRu)
                    pF[k] <- runif(1)
                    nbngbrs[k] <- runif(1, nbngbrsl, nbngbrsu)

                    S[, found_ind & (ftrial < S[1, ])] <- c(ftrial, trial)
                    if (sum(found_ind) > 1)
                        S <- unique(S, MARGIN = 2)
                } else if (ncol(S) < archive_size)
                    S <- cbind(S, c(ftrial, trial))
            }
        })
    } else {
        expression({
            if (ftrial < best_fpop ||
                isTRUE(all.equal(best_fpop, ftrial, tolerance = crit))) {
                if (is.null(S))
                    S <- as.matrix(c(ftrial, trial))
                if (ftrial < best_fpop)
                    best_fpop <- ftrial
                found_ind <- sqrt(colSums(
                    (trial - S[-1, , drop = FALSE])^2
                )) <= R
                if (any(found_ind)) {
                    S[, found_ind & (ftrial < S[1, ])] <- c(ftrial, trial)
                    if (sum(found_ind) > 1)
                        S <- unique(S, MARGIN = 2)
                } else if (ncol(S) < archive_size)
                    S <- cbind(S, c(ftrial, trial))
            }
        })
    }

    check_archive_constr <- if (reinit_if_solu_in_arch) {
        expression({
            if (all( htrial <= 0 ) && (ftrial < best_fpop ||
                isTRUE(all.equal(best_fpop, ftrial, tolerance = crit)))) {
                if (ftrial < best_fpop) {
                    if (is.null(S))
                        S <- as.matrix(c(ftrial, trial, htrial))
                    best_fpop <- ftrial
                }
                found_ind <- sqrt(colSums(
                    (trial - S[x_ind_in_S, , drop = FALSE])^2
                )) <= R
                if (any(found_ind)) {
                    # Re-initialize nearest neighbor of the trial vector
                    pop[, k] <- runif(d, lower, upper)
                    fpop[k] <- fn1(pop[, k])
                    hpop[, k] <- constr1(pop[, k])
                    F[, k] <- if (use_jitter)
                        runif(1, Fl, Fu) * (1 + jitter_factor*runif(d, -0.5, 0.5))
                    else runif(1, Fl, Fu)
                    CR[k] <- runif(1, CRl, CRu)
                    pF[k] <- runif(1)
                    nbngbrs[k] <- runif(1, nbngbrsl, nbngbrsu)
                    TAVpop[k] <- sum(pmax(hpop[, k], 0))

                    S[, found_ind & (ftrial < S[1, ])] <- c(ftrial, trial, htrial)
                    if (sum(found_ind) > 1)
                        S <- unique(S, MARGIN = 2)
                } else if (ncol(S) < archive_size)
                    S <- cbind(S, c(ftrial, trial, htrial))
            }
        })
    } else {
        expression({
            if (all( htrial <= 0 ) && (ftrial < best_fpop ||
                isTRUE(all.equal(best_fpop, ftrial, tolerance = crit)))) {
                if (ftrial < best_fpop) {
                    if (is.null(S))
                        S <- as.matrix(c(ftrial, trial, htrial))
                    best_fpop <- ftrial
                }
                found_ind <- sqrt(colSums(
                    (trial - S[x_ind_in_S, , drop = FALSE])^2
                )) <= R
                if (any(found_ind)) {
                    S[, found_ind & (ftrial < S[1, ])] <- c(ftrial, trial, htrial)
                    if (sum(found_ind) > 1)
                        S <- unique(S, MARGIN = 2)
                } else if (ncol(S) < archive_size)
                    S <- cbind(S, c(ftrial, trial, htrial))
            }
        })
    }

    identification_radius <- if (is.null(niche_radius)) {
        expression({
            dist <- vapply(
                pop_index,
                function(i) min(sqrt(colSums(
                    (pop[, i] - pop[, -i, drop = FALSE])^2
                ))),
                0
            )
            R <- min(R, mean(dist))
        })
    } else expression()

    update_pop <- if (is.null(constr)) {
        expression({
            pop <- pop_next
            fpop <- fpop_next
            F <- F_next
            CR <- CR_next
            pF <- pF_next
            nbngbrs <- nbngbrs_next
        })
    } else {
        expression({
            pop <- pop_next
            fpop <- fpop_next
            hpop <- hpop_next
            F <- F_next
            CR <- CR_next
            pF <- pF_next
            nbngbrs <- nbngbrs_next
            TAVpop <- TAVpop_next
        })
    }

    child <- if (is.null(constr)) { # Evaluate/select
        expression({
            ftrial <- fn1(trial)
            if (ftrial <= fpop[k]) {
                pop_next[, k] <- trial
                fpop_next[k] <- ftrial
                F_next[, k] <- Ftrial
                CR_next[k] <- CRtrial
                pF_next[k] <- pFtrial
                nbngbrs_next[k] <- nbngbrstrial
                eval(check_archive)
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
                if (TAVtrial <= TAVpop[k]) { # trial and target are both
                    pop_next[, k] <- trial   # unfeasible, the one with smaller
                    hpop_next[, k] <- htrial # constraint violation is chosen
                    F_next[, k] <- Ftrial    # or trial vector when both are
                    CR_next[k] <- CRtrial    # solutions of equal quality
                    pF_next[k] <- pFtrial
                    nbngbrs_next[k] <- nbngbrstrial
                    TAVpop_next[k] <- TAVtrial
                }
            } else if (TAVpop[k] > mu) { # trial is feasible and target is not
                pop_next[, k] <- trial
                fpop_next[k] <- fn1(trial)
                hpop_next[, k] <- htrial
                F_next[, k] <- Ftrial
                CR_next[k] <- CRtrial
                pF_next[k] <- pFtrial
                nbngbrs_next[k] <- nbngbrstrial
                TAVpop_next[k] <- TAVtrial
            } else {                       # between two feasible solutions, the
                ftrial <- fn1(trial)       # one with better objective function
                if (ftrial <= fpop[k]) {   # value is chosen
                    pop_next[, k] <- trial # or trial vector when both are
                    fpop_next[k] <- ftrial # solutions of equal quality
                    hpop_next[, k] <- htrial
                    F_next[, k] <- Ftrial
                    CR_next[k] <- CRtrial
                    pF_next[k] <- pFtrial
                    nbngbrs_next[k] <- nbngbrstrial
                    TAVpop_next[k] <- TAVtrial
                    eval(check_archive_constr)
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
                if (TAVtrial <= TAVpop[k]) { # trial and target both unfeasible
                    pop_next[, k] <- trial
                    hpop_next[, k] <- htrial
                    F_next[, k] <- Ftrial
                    CR_next[k] <- CRtrial
                    pF_next[k] <- pFtrial
                    nbngbrs_next[k] <- nbngbrstrial
                    TAVpop_next[k] <- TAVtrial
                }
            } else if (TAVpop[i] > mu) { # trial is feasible and target is not
                pop_next[, k] <- trial
                fpop_next[k] <- fn1(trial)
                hpop_next[, k] <- htrial
                F_next[, k] <- Ftrial
                CR_next[k] <- CRtrial
                pF_next[k] <- pFtrial
                nbngbrs_next[k] <- nbngbrstrial
                TAVpop_next[k] <- TAVtrial
                FF <- sum(TAVpop <= mu)/NP
                mu <- mu*(1 - FF/NP)
            } else {                       # two feasible solutions
                ftrial <- fn1(trial)
                if (ftrial <= fpop[k]) {
                    pop_next[, k] <- trial
                    fpop_next[k] <- ftrial
                    hpop_next[, k] <- htrial
                    F_next[, k] <- Ftrial
                    CR_next[k] <- CRtrial
                    pF_next[k] <- pFtrial
                    nbngbrs_next[k] <- nbngbrstrial
                    TAVpop_next[k] <- TAVtrial
                    eval(check_archive_constr)
                    FF <- sum(TAVpop <= mu)/NP
                    mu <- mu*(1 - FF/NP)
                }
            }
        })
    }

    fn1 <- function(par) fn(par, ...)

    if (!is.null(constr))
        constr1 <- if (meq > 0) {
            equalIndex <- 1:meq
            function(par) {
                h <- constr(par, ...)
                h[equalIndex] <- abs(h[equalIndex]) - eps
                h
            }
        } else function(par) constr(par, ...)

    use_jitter <- !is.null(jitter_factor)

    # Initialization
    d <- length(lower)
    pop <- matrix(runif(NP*d, lower, upper), nrow = d)
    if (!is.null(add_to_init_pop)) {
        pop <- unname(cbind(pop, add_to_init_pop))
        NP <- ncol(pop)
    }
    stopifnot(NP >= 4,
              length(nbngbrsl) == 1, is.numeric(nbngbrsl), nbngbrsl >= 3,
              length(nbngbrsu) == 1, is.numeric(nbngbrsu), nbngbrsu <= NP,
              nbngbrsl <= nbngbrsu)
    F <- if (use_jitter)
        (1 + jitter_factor*runif(d, -0.5, 0.5)) %o% runif(NP, Fl, Fu)
    else matrix(runif(NP, Fl, Fu), nrow = 1)
    CR <- runif(NP, CRl, CRu)
    pF <- runif(NP)
    nbngbrs <- runif(NP, nbngbrsl, nbngbrsu)
    fpop <- apply(pop, 2, fn1)
    stopifnot(is.vector(fpop), !anyNA(fpop), !is.nan(fpop))
    pop_next <- pop
    F_next <- F
    CR_next <- CR
    pF_next <- pF
    nbngbrs_next <- nbngbrs
    fpop_next <- fpop
    if (!is.null(constr)) {
        hpop <- apply(pop, 2, constr1)
        stopifnot(is.matrix(hpop) || is.vector(hpop),
                  !anyNA(hpop), !is.nan(hpop))
        if (is.vector(hpop)) dim(hpop) <- c(1, length(hpop))
        TAVpop <- apply( hpop, 2, function(x) sum(pmax(x, 0)) )
        mu <- median(TAVpop)
        hpop_next <- hpop
        TAVpop_next <- TAVpop
    }
    S <- NULL
    R <- if (is.null(niche_radius)) Inf else niche_radius
    best_fpop <- if (!is.null(constr)) Inf else min(fpop)
    x_ind_in_S <- 2:(d+1)

    pop_index <- 1:NP
    iteration <- 0

    while (iteration < maxiter) { # Generation loop
        iteration <- iteration + 1

        eval(identification_radius)

        for (i in pop_index) { # Start loop through population
            # Equalize the mean lifetime of all vectors
            # Price, KV, Storn, RM, and Lampinen, JA (2005)
            # Differential Evolution: A Practical Approach to
            # Global Optimization. Springer, p 284
            i <- ((iteration + i) %% NP) + 1

            # Fi update
            # Combine jitter with dither
            # Storn, Rainer (2008).
            # Differential evolution research - trends and open questions.
            # In: U. K. Chakraborty (Ed.), Advances in Differential Evolution,
            # SCI 143, Springer-Verlag, pp 11-12
            Ftrial <- if (runif(1) <= tau_F) {
                if (use_jitter)
                    runif(1, Fl, Fu) * (1 + jitter_factor*runif(d, -0.5, 0.5))
                else runif(1, Fl, Fu)
            } else F[, i]
            # CRi update
            CRtrial <- if (runif(1) <= tau_CR)
                runif(1, CRl, CRu)
            else CR[i]
            # pFi update
            pFtrial <- if (runif(1) <= tau_pF)
                runif(1)
            else pF[i]
            # nbngbrsi update
            nbngbrstrial <- if (runif(1) <= tau_nbngbrs)
                runif(1, nbngbrsl, nbngbrsu)
            else nbngbrs[i]

            # DE/rand/1/either-or/bin
            X_i <- pop[, i]
            # Select smallest distance members to the target vector
            subpop_ind <- order(
                sqrt(colSums((X_i - pop[, -i, drop = FALSE])^2))
            )[1:nbngbrstrial]
            # Randomly pick 3 vectors from the subpopulation
            r <- sample((pop_index[-i])[subpop_ind], 3)
            X_base <- pop[, r[1L]]
            X_r1 <- pop[, r[2L]]
            X_r2 <- pop[, r[3L]]

            trial <- handle_bounds(perform_reproduction(), X_base)
            # Identify the most similar individual of the trial vector
            k <- which.min( sqrt(colSums((trial - pop)^2)) )

            eval(child)
        }

        eval(update_pop)

        if (trace && (iteration %% triter == 0)) {
            x_best_in_pop <- which_best(fpop)
            x_best_in_S <- which.min(S[1, ])
            cat(iteration, ":", "<", R, ">",
                "population>>",
                "(", fpop[x_best_in_pop], ")", pop[, x_best_in_pop],
                if (!is.null(constr))
                    paste("{", which(hpop[, x_best_in_pop] > 0), "}"),
                "archive>>", "[", ncol(S), "]",
                "(", S[1, x_best_in_S], ")", S[x_ind_in_S, x_best_in_S],
                fill = TRUE)
        }
    }

    res <- list(iter = iteration)
    if (!is.null(S)) {
        ord <- order(S[1, ])
        res$solution_arch <- unname(S[x_ind_in_S, ord, drop = FALSE])
        res$objective_arch <- unname(S[1, ord])
        if (!is.null(constr))
            res$constr_value_arch <- unname(S[-(1:(d+1)), ord, drop = FALSE])
    }
    if (!is.null(constr)) {
        ord <- order(apply(hpop > 0, 2, any), fpop)
        res$constr_value_pop <- hpop[, ord, drop = FALSE]
    } else ord <- order(fpop)
    res$solution_pop <- pop[, ord, drop = FALSE]
    res$objective_pop <- fpop[ord]
    res
}
