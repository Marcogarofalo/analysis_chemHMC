#define functions_analysis_chemHMC_C
#include "functions_analysis_chemHMC.hpp"

double lhs_function_analysis_chemHMC_eg(int j, double**** in, int t, struct fit_type fit_info){
    double r=in[j][0][t][0];
    return r;
}
