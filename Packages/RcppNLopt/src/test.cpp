#include <Rcpp.h>
#include <nlopt.hpp>

using namespace Rcpp;


double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
  if (!grad.empty()) {
    grad[0] = 0.0;
    grad[1] = 0.5 / sqrt(x[1]);
  }
  return sqrt(x[1]);
}

typedef struct {
  double a, b;
} my_constraint_data;


double myconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data) {
  my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
  double a = d->a, b = d->b;
  if (!grad.empty()) {
    grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
    grad[1] = -1.0;
  }
  return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}



// [[Rcpp::export(rng=false)]]
NumericVector runTest() {

  // Create an instance of an optimization algorithm
  nlopt::opt opt(nlopt::LD_MMA, 2);

  // Set lower boundaries
  std::vector<double> lb{-HUGE_VAL, 0};
  //lb[0] = -HUGE_VAL; lb[1] = 0;
  opt.set_lower_bounds(lb);

  // Assign objective function
  opt.set_min_objective(myfunc, NULL);

  // Add inequality constraints
  my_constraint_data data[2] = { {2,0}, {-1,1} };
  opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
  opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);

  // Set convergence criterion
  opt.set_xtol_rel(1e-4);

  // Set x0
  std::vector<double> x{1.234, 5.678};

  // Run optimization. Solution is is passed back by reference.
  double minf;
  nlopt::result result = opt.optimize(x, minf);

  Rcout << result << '\n' << minf << '\n';

  return wrap(x);
}

