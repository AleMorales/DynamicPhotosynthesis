## WE have turn off closure checking

# context("closures")
#
# # Rosenbrock function and its gradient with original initial values
# # and parameter values (Rosenbrock, 1960, The Computer Journal)
# rosenbrock = function(par, a, b) {
#   x = par[[1]]
#   y = par[[2]]
#   result = (a - x)^2 + b*(y - x^2)^2
#   result
# }
#
# rosenbrock_grad = function(par, a, b) {
#   x = par[[1]]
#   y = par[[2]]
#   dfdx = -2*(a - x)  - 4*b*x*(y - x^2)
#   dfdy = 2*b*(y - x^2)
#   return(c(dfdx, dfdy))
# }
#
# rosenbrock2 = function(par) {
#   x = par[[1]]
#   y = par[[2]]
#   a = 1
#   b = 100
#   result = (a - x)^2 + b*(y - x^2)^2
#   result
# }
#
# rosenbrock_grad2 = function(par) {
#   x = par[[1]]
#   y = par[[2]]
#   a = 1
#   b = 100
#   dfdx = -2*(a - x)  - 4*b*x*(y - x^2)
#   dfdy = 2*b*(y - x^2)
#   return(c(dfdx, dfdy))
# }
#
#
# x0 = c(-1.2, 1)
# a = 1; b = 100
#
# test_that("wrong signatures", {
#   expect_error({nelder = nlopt(x0 = x0, algorithm = "nelder", ofun = "rosenbrock", a = a, b = b)})
#   expect_error({nelder = nlopt(x0 = x0, algorithm = "lbfgs", ofun = "rosenbrock_grad", a = a, b = b)})
# })
#
# test_that("ofun missing arguments", {
#   expect_error({nelder = nlopt(x0 = x0, algorithm = "nelder", ofun = rosenbrock, a = a)})
#   expect_error({nelder = nlopt(x0 = x0, algorithm = "nelder", ofun = rosenbrock)})
# })
#
# test_that("ofun extra arguments", {
#   expect_error({nelder = nlopt(x0 = x0, algorithm = "nelder", ofun = rosenbrock, a = a, b = b, c = 1)})
#   expect_error({nelder = nlopt(x0 = x0, algorithm = "nelder", ofun = rosenbrock2, a = a)})
# })
#
# test_that("grad missing arguments", {
#   expect_error({nelder = nlopt(x0 = x0, algorithm = "lbfgs", ofun = rosenbrock, grad = rosenbrock_grad, a = a)})
#   expect_error({nelder = nlopt(x0 = x0, algorithm = "lbfgs", ofun = rosenbrock, grad = rosenbrock_grad)})
# })
#
# test_that("grad extra arguments", {
#   expect_error({nelder = nlopt(x0 = x0, algorithm = "lbfgs", ofun = rosenbrock, grad = rosenbrock_grad, a = a, b = b, c = 1)})
#   expect_error({nelder = nlopt(x0 = x0, algorithm = "lbfgs", ofun = rosenbrock2, grad = rosenbrock_grad2, a = a)})
#
# })
