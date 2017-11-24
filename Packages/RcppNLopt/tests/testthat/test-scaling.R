context("scaling")


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



test_that("parscale", {
  expect_silent({nelder = nlopt(x0 = x0, algorithm = "nelder", ofun = rosenbrock, a = a, b = b,
                                options = list(parscale = abs(x0)))})
  expect_true(nelder$status > 0)
  expect_equal(nelder$par, c(a, a^2))
  expect_equal(nelder$min, 0)
})

test_that("fnscale", {
  expect_silent({nelder = nlopt(x0 = x0, algorithm = "nelder", ofun = function(x, a, b) -rosenbrock(x, a, b),
                                a = a, b = b, options = list(fnscale = -1))})
  expect_true(nelder$status > 0)
  expect_equal(nelder$par, c(a, a^2))
  expect_equal(nelder$min, 0)
})
