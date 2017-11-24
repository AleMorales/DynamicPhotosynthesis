context("bayfit")

# Log posterior. Estimate the mean of a normal distribution
lp = function(pars) {
  sum(dnorm(c(-10, 34), pars, 1, log = T)) + sum(dnorm(pars, 0, 20, log = T))
}

test_that("bayfit works", {
  fit = bayfit(lp = lp, start = c(-2, 10))
})
