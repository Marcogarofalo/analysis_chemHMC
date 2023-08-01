#define functions_analysis_chemHMC_C
#include "functions_analysis_chemHMC.hpp"

double rhs_const(int n, int Nvar, double *x, int Npar, double *P) {
  return P[0];
}

double lhs_radial_distr_function(int n, int e, int j, data_all gjack,
                                 struct fit_type fit_info) {
  double r;
  r = gjack.en[e].jack[fit_info.corr_id[0]][j];

  return r;
}
