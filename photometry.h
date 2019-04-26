/*
*Just some conviniece methods for converting between f_lambda and f_nu
*/

double C = 2.99792458e18;


double get_f_nu(double f_lam, double lambda) {
    return lambda*lambda*f_lam/C;
}

double get_nu(double lambda){
    return C/lambda;
}

double flux_den(double f_nu, double filter, double nu){
    return f_nu * filter/nu;
}

double filter_norm(double filter, double nu) {
    return filter/nu;
}
