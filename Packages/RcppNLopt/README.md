# RcppNLopt

The goal of RcppNLopt is to provide a high-level wrapper to the optimization library [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt). It supports all algorithms from NLopt (except Augmented Lagrangians and Multi-Level Single Linkage algorithms), for objective functions and constraints specified in R or C++ (using [Rcpp](http://www.rcpp.org)). In addition, it provides facilities for calculating  gradients of the objective function using differentiation methods.


## Calculation of gradients

Several optimization algorithms within NLopt require additional functions to calculate the gradient of the objective function and the functions implementing constraints on parameter space. If the user cannot provide such functions, RcppNLopt provides several methods to calculate or approximate such gradients. The methods available are:

* [Finite differences](https://en.wikipedia.org/wiki/Finite_difference): The gradients may be approximated numerically using forward, backward and central differences. This approach will work with any objective and constraint functions but but the numerical error in the approximation may not be sufficiently low for certain problems.

* [Symbolic differentiaton](https://en.wikipedia.org/wiki/Symbolic_computation): RcppNLopt will make use of the package [Deriv](https://cran.r-project.org/web/packages/Deriv) to calculate the gradient of the objective and constraint functions. This approach is limited as not all functions may be recognized by Deriv, but differentiation rules can be extended (see `?Deriv` for details). This approach generates the exact gradients and relies on source code transformation (i.e. new functions will be created to calculate the gradients).

* [Complex-step derivative approximation](https://en.wikipedia.org/wiki/Numerical_differentiation#Complex_variable_methods): This approach is similar to finite differences but relies on evaluating the objective and constraint functions with complex numbers as inputs. Thus, it is limited to functions that can be evaluated in the domain of complex numbers. It offers higher numerical precision than finite differences. Some issues include functions like `abs` which returns the modulus (the canonical interpretation, but useless in this case), and logical operators that require the use of `Re()` as they only operate on real-valued numbers (fortunately you can use `Re()` on a numeric variable and it has no effect, so the code will be compatible across real and complex  numbers).

* [Dual numbers](https://en.wikipedia.org/wiki/Automatic_differentiation#Automatic_differentiation_using_dual_numbers): This approach allows calculating the exact gradients for the objective and constraint functions and relies on evaluatng the objective and constraint functions with dual numbers as inputs. Thus, it is limited to functions that can be evaluated in the domain of dual numbers. This approach relies on the package [fad]().

The approach to implement these methods is for the user to call the functions `finite_grad`, `symbolic_grad`, `complex_grad` or `dual_grad` on the function for which the gradient is required and assign the result to the corresponding argument in the call to `nlopt`.
 
 
 
