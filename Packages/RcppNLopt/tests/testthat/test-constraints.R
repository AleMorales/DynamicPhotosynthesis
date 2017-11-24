context("constraints")


# Inequality constraint -------------------------------------------------------------------------------------------

# objective function
ofun = function(x, a, b) sqrt(x[2])

# gradient of objective function
ofun_grad = function(x, a, b) c(0, .5/sqrt(x[2]))

# constraint function
ineq = function(x, a, b) (a*x[1] + b)^3 - x[2]

# jacobian of constraint (a column for each contraint, number of rows = number of states)
ineq_grad = function( x, a, b )
  cbind(c(3*a[1]*(a[1]*x[1] + b[1])^2, -1.0),
        c(3*a[2]*(a[2]*x[1] + b[2])^2, -1.0))

# Data and initial values
a = c(2,-1)
b = c(0, 1)
x0=c(1,5)

# Solution
# minf = 0.54433104762009
# xopt = 0.3333333 0.2962963

test_that("inequality constraints work", {
  expect_silent({sqp = nlopt(x0 = x0, algorithm = "sqp",
                                lb = c(-Inf, 0), ub = c(Inf, Inf),
                                ofun = ofun, grad = ofun_grad,
                                ineq = ineq, ineq_grad = ineq_grad,
                                a = a, b = b)})
})


# Equality constraint -------------------------------------------------------------------------------------------

# objective function
ofun = function(x, params) 1

# gradient of objective function
ofun_grad = function(x, params) 0

# constraint function
eq = function(x, params) params[1]*x^2 + params[2]*x + params[3]

# jacobian of constraint (a column for each contraint, number of rows = number of states)
eq_grad = function(x, params)
  2*params[1]*x + params[2]

# Data and initial values
params = c(1, 1, -1)
x0 = -5

# Solution
# minf = 1
# xopt = -1.618033

test_that("inequality constraints work", {
  expect_silent({sqp = nlopt(x0 = x0, algorithm = "sqp",
                             ofun = ofun, grad = ofun_grad,
                             eq = eq, eq_grad = eq_grad,
                             params = params)})
})
