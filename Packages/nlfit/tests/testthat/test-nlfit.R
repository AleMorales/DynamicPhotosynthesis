context("nlfit")
library(nlfit)


dat = 1:100
Km = 25
Vmax = 50
y = rnorm(length(dat), Vmax*dat/(Km + dat), 5)
y[5] = NA

test_that("nlfit works", {
  fit = nlfit(model = function(pars) {pars[1]*dat/(pars[2] + dat)}, obs = y, start = c(V = 60, K = 15),
              algorithm = "subplex")
  fit2 = nlfit(model = y~V*x/(K + x), start = c(V = 60, K = 15),algorithm = "subplex",
                data = data.frame(dat = dat, y = y), lower = c(0,0))
})


prof_fit = profile(fit, lower = c(0,0), upper = c(K = 100), lower_delta = 0.005)

library(doParallel)
cl = makeCluster(2)
registerDoParallel(cl)
clusterExport(cl, "dat")
prof_fit = profile(fit, lower = c(0,0), upper = c(K = 100), lower_delta = 0.01, cl = cl)
stopCluster(cl)
plot(prof_fit)

confint(fit)
confint(fit, type = "hessian", hessian = "Richardson")
confint(fit, type = "hessian", hessian = "complex")
confint(fit, type = "hessian", hessian = "symbolic")
confint(prof_fit)
slice_fit = slice(fit)
plot(slice_fit)
confint(slice_fit)

prof_fit2 = profile(fit2)
plot(prof_fit2)
confint(fit2)
confint(fit2, type = "hessian", hessian = "Richardson")
confint(fit2, type = "hessian", hessian = "complex")
confint(fit2, type = "hessian", hessian = "symbolic")
confint(prof_fit2)

slice_fit2 = slice(fit2)
plot(slice_fit2)
confint(slice_fit2)
