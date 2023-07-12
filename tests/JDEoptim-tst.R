require("DEoptimR")

c.time <- function(...) cat('Time elapsed: ', ..., '\n')
S.time <- function(expr) c.time(system.time(expr))
source(system.file("xtraR/opt-test-funs.R", package = "DEoptimR"))
## sf1(), swf() + g11, RND, HEND, and alkylation list of $obj and $con testing functions
(doExtras <- DEoptimR:::doExtras())

set.seed(2345)
# Bound-constrained test problems ----------------------------------------------
S.time(sf1. <- JDEoptim(c(-100, -100), c(100, 100), sf1,
                        NP = 50, tol = 1e-7, maxiter = 800))
S.time(swf. <- JDEoptim(rep(-500, 10), rep(500, 10), swf,
                        tol = 1e-7))
# Only equality constraints ----------------------------------------------------
S.time(g11. <- JDEoptim(-c(1, 1), c(1, 1),
                        fn = g11$obj, constr = g11$con, meq = g11$eq, eps = 1e-3,
                        tol = 1e-7))
# Only inequality constraints --------------------------------------------------
S.time(RND. <- JDEoptim(c(1e-5, 1e-5), c(16, 16), RND$obj, RND$con,
                        NP = 40, tol = 1e-7))
if (doExtras) {
    S.time(HEND. <-
           JDEoptim(c(  100,  1000,  1000 ,  10,   10),
                    c(10000, 10000, 10000, 1000, 1000),
                    fn = HEND$obj, constr = HEND$con,
                    tol = 1e-4, trace = TRUE))
    S.time(alkylation. <-
           JDEoptim(c(1500,   1, 3000, 85, 90,  3, 145),
                    c(2000, 120, 3500, 93, 95, 12, 162),
                    fn = alkylation$obj, constr = alkylation$con,
                    tol = 0.1, trace = TRUE))
}

# Expected optimal values ------------------------------------------------------
bare.p.v <- function(r) unlist(unname(r[c("par", "value")]))
stopifnot(
    all.equal( bare.p.v(sf1.), c(0, 0, 0), tolerance = 1e-4 ),
    all.equal( bare.p.v(swf.), c(rep(420.97, 10), -418.9829*10),
               tolerance = 1e-4 ),
    all.equal( bare.p.v(RND.), c(3.036504, 5.096052, -0.388812),
               tolerance = 1e-2 ),
    all.equal( unname(c(abs(g11.$par[1]), g11.$par[2], g11.$value)),
               c(1/sqrt(2), 0.5, 0.75),
               tolerance = 1e-2 )
)
if (doExtras) {
    stopifnot(
        all.equal( bare.p.v(HEND.),
                   c(579.19, 1360.13, 5109.92, 182.01, 295.60, 7049.25),
                   tolerance = 1e-3 ),
        all.equal( bare.p.v(alkylation.),
                   c(1698.256922, 54.274463, 3031.357313, 90.190233,
                     95.0, 10.504119, 153.535355, -1766.36),
                   tolerance = 1e-2 )
    )
}

c.time(proc.time())
