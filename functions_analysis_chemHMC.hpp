#ifndef functions_analysis_chemHMC_H
#define functions_analysis_chemHMC_H
#include "non_linear_fit.hpp"

double rhs_const(int n, int Nvar, double* x, int Npar, double* P);
double lhs_radial_distr_function(int n, int e, int j, data_all gjack, struct fit_type fit_info);
#endif
