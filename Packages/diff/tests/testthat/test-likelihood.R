context("likelihood")

ofun = function(pars, data) {
  mu = exp(-pars[1]*data[,1])
  sigma = pars[2]
  -sum(dnorm(x = data[,1], mean = mu, sd = sigma, log = TRUE))
}
data = cbind(1:10, exp(-0.3*(1:10)))
sol = c(163.87047288175, -322.54497855496)

test_that("Gradient of NLL works", {
  expect_equal(sol, finite_grad(ofun, "Richardson")(c(0.1, 1), data = data))
  expect_equal(sol, finite_grad(ofun, "simple", h = 1e-8)(c(0.1, 1), data = data), tolerance = 1e-6)
  expect_equal(sol, finite_grad(ofun, "central", h = 1e-8)(c(0.1, 1), data = data))
  expect_equal(sol, dual_grad(ofun)(c(0.1, 1), data = data))
  expect_equal(sol, unname(symbolic_grad(ofun, c(pars = "1", pars = "2"))(c(0.1, 1), data = data)))
})


test_that("Gradients of NLL does not work", {
  expect_error(complex_jac(ofun)(c(0.1, 1), data = data), regexp = "function does not accept complex argument")
})

sol = rbind(c(-1154.17896935760, -327.74094577106), c(-327.74094577106, 987.63493567122))

test_that("Hessian of NLL works", {
  expect_equal(sol, finite_hess(ofun)(c(0.1, 1), data = data))
  expect_equal(sol, unname(symbolic_hess(ofun, c(pars = "1", pars = "2"))(c(0.1, 1), data = data)))
})
