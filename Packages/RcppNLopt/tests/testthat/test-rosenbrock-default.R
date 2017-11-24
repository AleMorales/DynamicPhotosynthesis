context("rosenbrock-defaults")


# Rosenbrock function and its gradient with original initial values
# and parameter values (Rosenbrock, 1960, The Computer Journal)
rosenbrock = function(par, a, b) {
  x = par[[1]]
  y = par[[2]]
  result = (a - x)^2 + b*(y - x^2)^2
  result
}

rosenbrock_grad = function(par, a, b) {
  x = par[[1]]
  y = par[[2]]
  dfdx = -2*(a - x)  - 4*b*x*(y - x^2)
  dfdy = 2*b*(y - x^2)
  return(c(dfdx, dfdy))
}


x0 = c(-1.2, 1)
a = 1; b = 100


# For each algorith test:
#  1. That it runs without producing any messages, errors or warnings
#  2. That the solution is found
# Issues with the tests:
# 1. For some algorithms the tolerance on expect_equal() must be increased above sqrt(epsilon)
# 2. Some algorithms need maxiter > 1000 to produce reasonable results and to relax the tolerance
#       ESCH is particularly bad
# 3. Global algorithms need boundaries (reasonable values are used)


# Local, derivative-free --------------------------------------------------

test_that("rosenbrock-nelder", {
  expect_silent({nelder = nlopt(x0 = x0, algorithm = "nelder", ofun = rosenbrock, a = a, b = b)})
  expect_true(nelder$status > 0)
  expect_equal(nelder$par, c(a, a^2))
  expect_equal(nelder$min, 0)
})

test_that("rosenbrock-subplex", {
  expect_silent({subplex = nlopt(x0 = x0, algorithm = "subplex", ofun = rosenbrock, a = a, b = b)})
  expect_true(subplex$status > 0)
  expect_equal(subplex$par, c(a, a^2), tolerance = 1.e-7)
  expect_equal(subplex$min, 0)
})

test_that("rosenbrock-newuoa", {
  expect_silent({newuoa = nlopt(x0 = x0, algorithm = "newuoa", ofun = rosenbrock, a = a, b = b)})
  expect_true(newuoa$status > 0)
  expect_equal(newuoa$par, c(a, a^2), tolerance = 1.5e-5)
  expect_equal(newuoa$min, 0)
})

test_that("rosenbrock-cobyla", {
  expect_silent({cobyla = nlopt(x0 = x0, algorithm = "cobyla", ofun = rosenbrock, a = a, b = b,
                                options = list(maxiter = Inf))})
  expect_true(cobyla$status > 0)
  expect_equal(cobyla$par, c(a, a^2), tolerance = 1e-4)
  expect_equal(cobyla$min, 0)
})

test_that("rosenbrock-bobyqa", {
  expect_silent({bobyqa = nlopt(x0 = x0, algorithm = "bobyqa", ofun = rosenbrock, a = a, b = b)})
  expect_true(bobyqa$status > 0)
  expect_equal(bobyqa$par, c(a, a^2))
  expect_equal(bobyqa$min, 0)
})

test_that("rosenbrock-praxis", {
  expect_silent({praxis = nlopt(x0 = x0, algorithm = "praxis", ofun = rosenbrock, a = a, b = b)})
  expect_true(praxis$status > 0)
  expect_equal(praxis$par, c(a, a^2))
  expect_equal(praxis$min, 0)
})


# Local, derivative-based -------------------------------------------------

test_that("rosenbrock-lbfgs", {
  expect_silent({lbfgs = nlopt(x0 = x0, algorithm = "lbfgs", ofun = rosenbrock, a = a, b = b,
                               grad = rosenbrock_grad)})
  expect_true(lbfgs$status > 0)
  expect_equal(lbfgs$par, c(a, a^2))
  expect_equal(lbfgs$min, 0)
})

test_that("rosenbrock-sqp", {
  expect_silent({sqp = nlopt(x0 = x0, algorithm = "sqp", ofun = rosenbrock, a = a, b = b,
                               grad = rosenbrock_grad)})
  expect_true(sqp$status > 0)
  expect_equal(sqp$par, c(a, a^2))
  expect_equal(sqp$min, 0)
})

test_that("rosenbrock-mma", {
  expect_silent({mma = nlopt(x0 = x0, algorithm = "mma", ofun = rosenbrock, a = a, b = b,
                             grad = rosenbrock_grad, options = list(maxiter = Inf))})
  expect_true(mma$status > 0)
  expect_equal(mma$par, c(a, a^2), tolerance = 3e-5)
  expect_equal(mma$min, 0)
})

test_that("rosenbrock-var1", {
  expect_silent({var = nlopt(x0 = x0, algorithm = "var", ofun = rosenbrock, a = a, b = b,
                             grad = rosenbrock_grad)})
  expect_true(var$status > 0)
  expect_equal(var$par, c(a, a^2))
  expect_equal(var$min, 0)
})

test_that("rosenbrock-var2", {
  expect_silent({var = nlopt(x0 = x0, algorithm = "var", ofun = rosenbrock, a = a, b = b,
                             grad = rosenbrock_grad, options = list(rank = 2))})
  expect_true(var$status > 0)
  expect_equal(var$par, c(a, a^2))
  expect_equal(var$min, 0)
})

test_that("rosenbrock-tnewton", {
  expect_silent({tnewton = nlopt(x0 = x0, algorithm = "tnewton", ofun = rosenbrock, a = a, b = b,
                             grad = rosenbrock_grad)})
  expect_true(tnewton$status > 0)
  expect_equal(tnewton$par, c(a, a^2))
  expect_equal(tnewton$min, 0)
})

test_that("rosenbrock-tnewton_restart", {
  expect_silent({tnewton_restart = nlopt(x0 = x0, algorithm = "tnewton", ofun = rosenbrock, a = a, b = b,
                                 grad = rosenbrock_grad, options = list(restart = TRUE))})
  expect_true(tnewton_restart$status > 0)
  expect_equal(tnewton_restart$par, c(a, a^2))
  expect_equal(tnewton_restart$min, 0)
})

test_that("rosenbrock-tnewton_precondition", {
  expect_silent({tnewton_precondition = nlopt(x0 = x0, algorithm = "tnewton", ofun = rosenbrock, a = a, b = b,
                                         grad = rosenbrock_grad, options = list(precondition = TRUE))})
  expect_true(tnewton_precondition$status > 0)
  expect_equal(tnewton_precondition$par, c(a, a^2))
  expect_equal(tnewton_precondition$min, 0)
})

test_that("rosenbrock-tnewton_precondition_restart", {
  expect_silent({tnewton_precondition_restart = nlopt(x0 = x0, algorithm = "tnewton", ofun = rosenbrock, a = a, b = b,
                                              grad = rosenbrock_grad, options = list(precondition = TRUE, restart = TRUE))})
  expect_true(tnewton_precondition_restart$status > 0)
  expect_equal(tnewton_precondition_restart$par, c(a, a^2))
  expect_equal(tnewton_precondition_restart$min, 0)
})


# Global, derivative-free -------------------------------------------------

test_that("rosenbrock-direct", {
  expect_silent({direct = nlopt(x0 = x0, algorithm = "direct", ofun = rosenbrock, a = a, b = b,
                                lb = c(-10, -10), ub = c(10,10),
                                options = list(maxiter = 2e4, scaling = TRUE, local = FALSE))})
  expect_true(direct$status > 0)
  expect_equal(direct$par, c(a, a^2))
  expect_equal(direct$min, 0)
})

test_that("rosenbrock-direct-local", {
  expect_silent({direct = nlopt(x0 = x0, algorithm = "direct", ofun = rosenbrock, a = a, b = b,
                                lb = c(-10, -10), ub = c(10,10),
                                options = list(maxiter = 1e4, local = TRUE, scaling = TRUE))})
  expect_true(direct$status > 0)
  expect_equal(direct$par, c(a, a^2))
  expect_equal(direct$min, 0)
})

test_that("rosenbrock-direct-local-random", {
  expect_silent({direct = nlopt(x0 = x0, algorithm = "direct", ofun = rosenbrock, a = a, b = b,
                                lb = c(-10, -10), ub = c(10,10),
                                options = list(maxiter = 1e4, local = TRUE, randomize = TRUE, scaling = TRUE))})
  expect_true(direct$status > 0)
  expect_equal(direct$par, c(a, a^2))
  expect_equal(direct$min, 0)
})

test_that("rosenbrock-direct-noscal", {
  expect_silent({direct = nlopt(x0 = x0, algorithm = "direct", ofun = rosenbrock, a = a, b = b,
                                lb = c(-10, -10), ub = c(10,10),
                                options = list(maxiter = 1e4, scaling = FALSE, local = FALSE))})
  expect_true(direct$status > 0)
  expect_equal(direct$par, c(a, a^2))
  expect_equal(direct$min, 0)
})


test_that("rosenbrock-direct-noscal-local", {
  expect_silent({direct = nlopt(x0 = x0, algorithm = "direct", ofun = rosenbrock, a = a, b = b,
                                lb = c(-10, -10), ub = c(10,10),
                                options = list(maxiter = 1e4, scaling = FALSE, local = TRUE))})
  expect_true(direct$status > 0)
  expect_equal(direct$par, c(a, a^2))
  expect_equal(direct$min, 0)
})

test_that("rosenbrock-direct-noscal-local-random", {
  expect_silent({direct = nlopt(x0 = x0, algorithm = "direct", ofun = rosenbrock, a = a, b = b,
                                lb = c(-10, -10), ub = c(10,10),
                                options = list(maxiter = 1e4, scaling = FALSE, local = TRUE, randomize = TRUE))})
  expect_true(direct$status > 0)
  expect_equal(direct$par, c(a, a^2))
  expect_equal(direct$min, 0)
})

test_that("rosenbrock-crs", {
  expect_silent({crs = nlopt(x0 = x0, algorithm = "crs", ofun = rosenbrock, a = a, b = b,
                               lb = c(-10, -10), ub = c(10,10),
                               options = list(maxiter = 1e5, randseed = 0))})
  expect_true(crs$status > 0)
  expect_equal(crs$par, c(a, a^2), tolerance = 1e-7)
  expect_equal(crs$min, 0)
})

test_that("rosenbrock-isres", {
  expect_silent({isres = nlopt(x0 = x0, algorithm = "isres", ofun = rosenbrock, a = a, b = b,
                             lb = c(-10, -10), ub = c(10,10),
                             options = list(maxiter = 5e4, randseed = 0))})
  expect_true(isres$status > 0)
  expect_equal(isres$par, c(a, a^2), tolerance = 1.5e-7)
  expect_equal(isres$min, 0)
})


test_that("rosenbrock-esch", {
  expect_silent({esch = nlopt(x0 = x0, algorithm = "esch", ofun = rosenbrock, a = a, b = b,
                               lb = c(-10, -10), ub = c(10,10),
                               options = list(maxiter = 1e5, randseed = 0))})
  expect_true(esch$status > 0)
  expect_equal(esch$par, c(a, a^2), tolerance = 0.3)
  expect_equal(esch$min, 0, tolerance = 5e-3)
})


# Global, derivative-based ------------------------------------------------

test_that("rosenbrock-stogo", {
  expect_silent({stogo = nlopt(x0 = x0, algorithm = "stogo", ofun = rosenbrock, a = a, b = b,
                                lb = c(-10, -10), ub = c(10,10), grad = rosenbrock_grad,
                                options = list(maxiter = 1e4, randomize = FALSE))})
  expect_true(stogo$status > 0)
  expect_equal(stogo$par, c(a, a^2), tolerance = 1e-7)
  expect_equal(stogo$min, 0)
})

test_that("rosenbrock-stogo-random", {
  expect_silent({stogo = nlopt(x0 = x0, algorithm = "stogo", ofun = rosenbrock, a = a, b = b,
                               lb = c(-10, -10), ub = c(10,10), grad = rosenbrock_grad,
                               options = list(maxiter = 1e4, randomize = TRUE))})
  expect_true(stogo$status > 0)
  expect_equal(stogo$par, c(a, a^2), tolerance = 1e-7)
  expect_equal(stogo$min, 0)
})
