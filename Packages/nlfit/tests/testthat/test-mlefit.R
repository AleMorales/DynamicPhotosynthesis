context("mlefit")


x = 1:100
Km = 25
Vmax = 50
y = rnorm(length(x), Vmax*x/(Km + x), 5)
ofun = function(pars) {
  mod = pars[1]*x/(pars[2] + x)
  -sum(-log(sqrt(2*pi)) - log(pars[3]) - (y - mod)^2/(2*pars[3]^2))
}

test_that("mlefit works", {
  fit = mlefit(NLL = ofun, start = c(V = 60, K = 15, sigma = 1))
})

slice_fit = slice(fit)
profile_fit = profile(fit)

plot(profile_fit)
plot(slice_fit)

confint(profile_fit)
confint(slice_fit)

confint(fit, type = "hessian", hessian = "Richardson")
confint(fit, type = "hessian", hessian = "complex")
confint(fit, type = "profile")
confint(fit, type = "slice")
