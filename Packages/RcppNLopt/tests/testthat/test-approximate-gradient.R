context("approximate-gradient")


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


forward_rosenbrock_grad = finite_grad(ofun = rosenbrock, method = "forward")
central_rosenbrock_grad = finite_grad(ofun = rosenbrock, method = "central")
complex_rosenbrock_grad = complex_grad(ofun = rosenbrock)
symbolic_rosenbrock_grad = symbolic_grad(ofun = rosenbrock, c("x", "y"))
fad_rosenbrock_grad = dual_grad(ofun = rosenbrock)

x0 = c(-1.2, 1)
a = 1; b = 100




test_that("differentiation", {
  expect_equal(forward_rosenbrock_grad(x0, a, b), rosenbrock_grad(x0, a, b), tolerance = 5e-6)
  expect_equal(central_rosenbrock_grad(x0, a, b), rosenbrock_grad(x0, a, b), tolerance = 1e-8)
  expect_equal(complex_rosenbrock_grad(x0, a, b), rosenbrock_grad(x0, a, b), tolerance = 1e-8)
  expect_equal(unname(symbolic_rosenbrock_grad(x0, a, b)), rosenbrock_grad(x0, a, b), tolerance = 1e-8)
  expect_equal(as.numeric(fad_rosenbrock_grad(x0, a, b)), rosenbrock_grad(x0, a, b), tolerance = 1e-8)
})

