require(DEoptimR)

c.time <- function(...) cat('Time elapsed: ', ..., '\n')
S.time <- function(expr) c.time(system.time(expr))
(doExtras <- DEoptimR:::doExtras())

set.seed(2345)
# Bound-constrained test problems ----------------------------------------------

bl <- function(x) {
#   Becker and Lago problem
#
#   -10 <= x1, x2 <= 10
#   The function has four minima located at (+-5, +-5), all with f(x*) = 0.
#
#   Source:
#     Ali, M. Montaz, Khompatraporn, Charoenchai, and Zabinsky, Zelda B. (2005).
#     A numerical evaluation of several stochastic algorithms on selected
#     continuous global optimization test problems.
#     Journal of Global Optimization 31, 635-672.

    sum((abs(x) - 5)^2)
}

S.time(bl_ <- NCDEoptim(-c(10, 10), c(10, 10),
                        bl,
                        niche_radius = 5,
                        maxiter = 100))
# Only inequality constraints --------------------------------------------------

#   Function F1
#
#   f(x) = x^2
#   subject to:
#   g(x) = 1 - x^2 <= 0
#
#   -2 <= x <= 2
#   The two global optima are (x1*, x2*; f*) = (1, -1; 1).
#
#   Source:
#     Poole, Daniel J. and Allen, Christian B. (2019).
#     Constrained niching using differential evolution.
#     Swarm and Evolutionary Computation 44, 74-100.

S.time(F1_ <- NCDEoptim(-2, 2,
                        \(x) x^2,
                        \(x) 1 - x^2,
                        niche_radius = 1,
                        maxiter = 200))

# Expected optimal values ------------------------------------------------------

stopifnot(
  all.equal( as.vector(abs(bl_$solution_arch)), rep(5, 8), tolerance = 1e-3 ),
  all.equal( as.vector(abs(F1_$solution_arch)), c(1, 1), tolerance = 1e-2 )
)

c.time(proc.time())
