/*
 *
 * Dependencies and selective imports
 *
 */
#include <vector>
#include <string>
#include <stdexcept>
#include <Rcpp.h>
#include <nlopt.hpp>

using std::vector;
using std::string;
using std::unordered_map;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::Function;
using Rcpp::List;
using Rcpp::wrap;
using Rcpp::as;
using Rcpp::_;
using Rcpp::Nullable;
using Rcpp::Rcout;
/*
 *
 * Structure that stores R closures and interfaces to STL
 *
 */
struct Rdata {

  // R functions
  const Function& Rofun;
  const Nullable<Function>& Rgrad;
  const Nullable<Function>& Rineq;
  const Nullable<Function>& Rgrad_ineq;
  const Nullable<Function>& Req;
  const Nullable<Function>& Rgrad_eq;
  int nevals = 0;

  Rdata(const Function& iRofun,
        const Nullable<Function>& iRgrad,
        const Nullable<Function>& iRineq,
        const Nullable<Function>& iRgrad_ineq,
        const Nullable<Function>& iReq,
        const Nullable<Function>& iRgrad_eq) :
    Rofun{iRofun}, Rgrad{iRgrad},
    Rineq{iRineq}, Rgrad_ineq{iRgrad_ineq},
    Req{iReq}, Rgrad_eq{iRgrad_eq} {};

  // Interfaces to STL
  double ofun(const vector<double>& x) {
    try {
      NumericVector result =  Rofun(wrap(x));
      nevals += 1;
      return result[0];
    } catch(const std::exception& e) {
      Rcout << "Error inside the ofun function: " << e.what() << '\n'; // Needed because nlopt will not print the message
      throw e;
    }
  };

  void grad(const vector<double>& x, vector<double>& grad) {
    try {
      NumericVector result = as<Function>(Rgrad)(wrap(x));
      grad = as<vector<double>>(result);
    } catch(const std::exception& e) {
      Rcout << "Error inside the grad function: " << e.what() << '\n'; // Needed because nlopt will not print the message
      throw e;
    }
  };

  NumericVector ineq(const NumericVector& x) {
    try {
      NumericVector result = as<Function>(Rineq)(x);
      return result;
    } catch(const std::exception& e) {
      Rcout << "Error inside the ineq function: " << e.what() << '\n'; // Needed because nlopt will not print the message
      throw e;
    }
  };
  void grad_ineq(const NumericVector& x, double* grad) {
    try {
      NumericMatrix result = as<Function>(Rgrad_ineq)(x);
      unsigned n = result.nrow(); //nstates
      unsigned m = result.ncol(); //ncontraints
      for(unsigned i = 0; i < m; i++)
        for(unsigned j = 0; j < n; j++)
          grad[i*n + j] = result(i,j);
    } catch(const std::exception& e) {
      Rcout << "Error inside the grad_ineq function: " << e.what() << '\n'; // Needed because nlopt will not print the message
      throw e;
    }
  };

  NumericVector eq(const NumericVector& x) {
    try {
      NumericVector result = as<Function>(Req)(x);
      return result;
    } catch(const std::exception& e) {
      Rcout << "Error inside the eq function: " << e.what() << '\n'; // Needed because nlopt will not print the message
      throw e;
    }
  };
  void grad_eq(const NumericVector& x, double* grad) {
    try {
      NumericMatrix result = as<Function>(Rgrad_eq)(x);
      unsigned n = result.nrow(); //nstates
      unsigned m = result.ncol(); //ncontraints
      for(unsigned i = 0; i < m; i++)
        for(unsigned j = 0; j < n; j++)
          grad[i*n + j] = result(i,j);
    } catch(const std::exception& e) {
      Rcout << "Error inside the grad_eq function: " << e.what() << '\n'; // Needed because nlopt will not print the message
      throw e;
    }
  };

};



/*
 *
 * Wrappers between NLopt and closures defined in R
 *
 */

// Callable from nlopt: fun and grad_fun in one single function
double ofun(const vector<double>& x, vector<double>& grad, void* f_data) {

  Rdata* data = static_cast<Rdata*>(f_data);

  double result = data->ofun(x);

  if(!grad.empty()) data->grad(x, grad);

  return result;
}

// Callable from nlopt: inequality constraints
void ineq(unsigned m, double *result_, unsigned n, const double* x_, double* grad_, void* f_data){
  // Convert input form the callback
  Rdata* data = static_cast<Rdata*>(f_data);
  NumericVector x(n);
  for(unsigned i = 0; i < n; i++) x[i] = x_[i];

  // Store the results of the constraints in pointer
  NumericVector result = data->ineq(x);
  if(m != result.size()) throw std::runtime_error("The number of constraints and length of ineqatol do not coincide");
  for(unsigned i = 0; i < m; i++) result_[i] = result[i];

  if(grad_) data->grad_ineq(x, grad_);
}


// Callable from nlopt: equality constraints
void eq(unsigned m, double *result_, unsigned n, const double* x_, double* grad_, void* f_data){
  // Convert input form the callback
  Rdata* data = static_cast<Rdata*>(f_data);
  NumericVector x(n);
  for(unsigned i = 0; i < n; i++) x[i] = x_[i];

  // Store the results of the constraints in pointer
  NumericVector result = data->eq(x);
  if(m != result.size()) throw std::runtime_error("The number of constraints and length of eqatol do not coincide");
  for(unsigned i = 0; i < m; i++) result_[i] = result[i];

  if(grad_) data->grad_eq(x, grad_);
}

// Because I cannot pass an enum from R, map string to enum
const static unordered_map<string, nlopt::algorithm> nlopt_algorithms {
  {"GN_DIRECT", nlopt::GN_DIRECT},
  {"GN_DIRECT_L", nlopt::GN_DIRECT_L},
  {"GN_DIRECT_L_RAND", nlopt::GN_DIRECT_L_RAND},
  {"GN_DIRECT_NOSCAL", nlopt::GN_DIRECT_NOSCAL},
  {"GN_DIRECT_L_NOSCAL", nlopt::GN_DIRECT_L_NOSCAL},
  {"GN_DIRECT_L_RAND_NOSCAL", nlopt::GN_DIRECT_L_RAND_NOSCAL},
  {"GD_STOGO", nlopt::GD_STOGO},
  {"GD_STOGO_RAND", nlopt::GD_STOGO_RAND},
  {"LD_LBFGS", nlopt::LD_LBFGS},
  {"LN_PRAXIS", nlopt::LN_PRAXIS},
  {"LD_VAR1", nlopt::LD_VAR1},
  {"LD_VAR2", nlopt::LD_VAR2},
  {"LD_TNEWTON", nlopt::LD_TNEWTON},
  {"LD_TNEWTON_RESTART", nlopt::LD_TNEWTON_RESTART},
  {"LD_TNEWTON_PRECOND", nlopt::LD_TNEWTON_PRECOND},
  {"LD_TNEWTON_PRECOND_RESTART", nlopt::LD_TNEWTON_PRECOND_RESTART},
  {"GN_CRS2_LM", nlopt::GN_CRS2_LM},
  {"GN_MLSL", nlopt::GN_MLSL},
  {"GD_MLSL", nlopt::GD_MLSL},
  {"GN_MLSL_LDS", nlopt::GN_MLSL_LDS},
  {"GD_MLSL_LDS", nlopt::GD_MLSL_LDS},
  {"LD_MMA", nlopt::LD_MMA},
  {"LN_COBYLA", nlopt::LN_COBYLA},
  {"LN_NEWUOA_BOUND", nlopt::LN_NEWUOA_BOUND},
  {"LN_NELDERMEAD", nlopt::LN_NELDERMEAD},
  {"LN_SBPLX", nlopt::LN_SBPLX},
  {"LN_AUGLAG", nlopt::LN_AUGLAG},
  {"LD_AUGLAG", nlopt::LD_AUGLAG},
  {"LN_AUGLAG_EQ", nlopt::LN_AUGLAG_EQ},
  {"LD_AUGLAG_EQ", nlopt::LD_AUGLAG_EQ},
  {"LN_BOBYQA", nlopt::LN_BOBYQA},
  {"GN_ISRES", nlopt::GN_ISRES},
  {"AUGLAG", nlopt::AUGLAG},
  {"AUGLAG_EQ", nlopt::AUGLAG_EQ},
  {"G_MLSL", nlopt::G_MLSL},
  {"G_MLSL_LDS", nlopt::G_MLSL_LDS},
  {"LD_SLSQP", nlopt::LD_SLSQP},
  {"LD_CCSAQ", nlopt::LD_CCSAQ},
  {"GN_ESCH", nlopt::GN_ESCH}
};

// Translate flag to message
const static unordered_map<int, string> nlopt_messages {
  {-1, "Generic failure code."},
  {-2, "Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera)."},
  {-3,"Ran out of memory."},
  {-4, "Halted because roundoff errors limited progress."},
  {-5, "Forced termination: nlopt_force_stop() was called."},
  {1, "Generic success return value."},
  {2, "Optimization stopped because stopval was reached."},
  {3, "Optimization stopped because ftol_rel or ftol_abs was reached."},
  {4, "Optimization stopped because xtol_rel or xtol_abs (above) was reached."},
  {5, "Optimization stopped because maxeval (above) was reached."},
  {6, "Optimization stopped because maxtime (above) was reached."}
};



// Only Rofun must be passed to C++, the rest of functions are by default NULL
// If the problem has no constraints, this is detected from the NULLity of the argument
// It is the responsibility of the R code to pass gradient functions for algorithms that need them
// If a gradient is missing, Rcpp will raise the exception "cannot convert to function"
// [[Rcpp::export(rng=false)]]
List wrap_nlopt(const NumericVector& x0,
                const NumericVector& lb,
                const NumericVector& ub,
                const List& options,
                const Function& Rofun,
                const Nullable<Function>& Rgrad = R_NilValue,
                const Nullable<Function>& Rineq = R_NilValue,
                const Nullable<Function>& Rgrad_ineq = R_NilValue,
                const Nullable<Function>& Req = R_NilValue,
                const Nullable<Function>& Rgrad_eq = R_NilValue) {

  // Create object that stores the R functions
  Rdata Rfuns{Rofun, Rgrad, Rineq, Rgrad_ineq, Req, Rgrad_eq};

  // Create nlopt optimizer object
  // Convert from string to enum and check algorithm exists
  auto Ralg = as<string>(options["nlopt_algorithm"]);
  auto alg = nlopt_algorithms.find(Ralg);
  if(alg == nlopt_algorithms.end())
    throw std::runtime_error("Algorithm " +  Ralg + " does not exist");
  nlopt::opt optim(alg->second, x0.size());

  // Use ofun wrapper to interface with nlopt "void* callback" style
  optim.set_min_objective(ofun, &Rfuns);

  // Same "void* callback" style for constraints. Only add when not null.
  // Only support simple constraints with scalar absolute tolerance
  if(Rineq.isNotNull())
    optim.add_inequality_mconstraint(ineq, &Rfuns, as<vector<double>>(options["ineqatol"]));
  if(Req.isNotNull())
    optim.add_equality_mconstraint(eq, &Rfuns, as<vector<double>>(options["eqatol"]));

  // Set lower and upper bounds of the problem
  optim.set_lower_bounds(as<vector<double>>(lb));
  optim.set_upper_bounds(as<vector<double>>(ub));

  // Termination criteria (options should contain valid values for all the criteria)
  NumericVector xatol = options["xatol"];
  if(xatol.size() > 1) {
    optim.set_xtol_abs(as<vector<double>>(xatol));
  } else {
    optim.set_xtol_abs(xatol[0]);
  }
  optim.set_xtol_rel(as<double>(options["xrtol"]));
  optim.set_ftol_abs(as<double>(options["fatol"]));
  optim.set_ftol_rel(as<double>(options["frtol"]));
  optim.set_maxeval(as<double>(options["maxiter"]));
  optim.set_maxtime(as<double>(options["maxtime"]));

  // Initial step for df algorithms (dx should not be present if derivative-based)
  if(options.containsElementNamed("dx")) {
    NumericVector dx = options["dx"];
    if(dx.size() > 1) {
      optim.set_initial_step(as<vector<double>>(dx));
    } else if(dx.size() == 1){
      optim.set_initial_step(dx[0]);
    }
  }

  // Pass random seed to nlopt if needed otherwise use time as seed.
  // randseed is only define for algorithms it matters,
  // but we need to use empty vector to signal we don't want to use it
  if(options.containsElementNamed("randseed")) {
    Rcpp::IntegerVector randseed = as<Rcpp::IntegerVector>(options["randseed"]);
    if(randseed.size() > 0) {
      nlopt::srand(randseed[0]);
    } else {
      nlopt::srand_time();
    }
  } else {
    nlopt::srand_time();
  }

  // User-defined population size (optional)
  if(options.containsElementNamed("pop")) {
    optim.set_population(as<int>(options["pop"]));
  }

  // User-defined M (number of gradients remembered by truncated Newton scheme)
  if(options.containsElementNamed("M")) {
    optim.set_vector_storage(as<int>(options["M"]));
  }


  // Pass options for local optimizers
  if(options.containsElementNamed("local_opts")) {
    List local_opts = options["local_opts"];

    // Create optimization object for local optimizer
    auto Ralg = as<string>(options["nlopt_algorithm"]);
    auto alg = nlopt_algorithms.find(Ralg);
    if(alg == nlopt_algorithms.end())
      throw std::runtime_error("Algorithm " +  Ralg + " does not exist");
    nlopt::opt local_optim(alg->second, x0.size());

    // Assign all the options except the functions...
    NumericVector xatol = options["xatol"];
    if(xatol.size() > 1) {
      local_optim.set_xtol_abs(as<vector<double>>(xatol));
    } else {
      local_optim.set_xtol_abs(xatol[0]);
    }
    local_optim.set_xtol_rel(as<double>(options["xrtol"]));
    local_optim.set_ftol_abs(as<double>(options["fatol"]));
    local_optim.set_ftol_rel(as<double>(options["frtol"]));
    local_optim.set_maxeval(as<double>(options["maxiter"]));
    local_optim.set_maxtime(as<double>(options["maxtime"]));

    // Assign local options
    optim.set_local_optimizer(local_optim);

  }


  // Actual optimization (results returned by reference status flag as output).
  double fmin = 0;
  vector<double> xsol = as<vector<double>>(x0);
  auto oflag = static_cast<int>(optim.optimize(xsol, fmin));

  // Organize output into a list
  // Only add info that can only be gathered in C++ (info on inputs better managed in R)
  List solution = List::create(
    _["min"] = wrap(fmin),
    _["par"] = wrap(xsol),
    _["status"] = wrap(oflag),
    _["message"] = wrap(nlopt_messages.find(oflag)->second),
    _["nevals"] = wrap(Rfuns.nevals)
  );

  return solution;
}




