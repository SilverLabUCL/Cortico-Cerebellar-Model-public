#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <gsl/gsl_sf_erf.h>
#include <numeric>
#include <limits>

///////////////////////////////////
// Heaviside function
///////////////////////////////////
inline double heavi(double x) {
	return (x >= 0.0) ? 1.0 : 0.0;
}

///////////////////////////////////
// tanh(x) function
///////////////////////////////////
inline double phi_tanh(double x, double a){
   		return  std::tanh(x);
}

///////////////////////////////////
// linear function
///////////////////////////////////
inline double phi_linear(double x, double a){
   		return  x;
}

///////////////////////////////////
// linear rectifier function
///////////////////////////////////
inline double phi_relu(double x, double a){
   		return  std::max(x,0.0);
}

///////////////////////////////////
// log(1+exp(x)) function
///////////////////////////////////
inline double phi_softplus(double x, double a) {
    double z = a * x;
    if (z > 0) {
        return (z + std::log1p(std::exp(-z))) / a;
    } else {
        return std::log1p(std::exp(z)) / a;
    }
}

///////////////////////////////////
// exp(x)/(1+exp(x)) function
///////////////////////////////////
inline double inv_logit(double x, double a) {
    double z = a * x;
    const double LIMIT = 40.0;
	if (z > LIMIT) return x;
	if (z < -LIMIT) return std::exp(z) / a;
    if (z >= 0) {
        double ez = std::exp(-z);
        return 1.0 / (1.0 + ez);
    } else {
        double ez = std::exp(z);
        return ez / (1.0 + ez);
    }
}

///////////////////////////////////
// nearestSPD : computes the the nearest (in Frobenius norm) Symmetric Positive Definite matrix to RHO
///////////////////////////////////
void nearestSPD(arma::mat& RHO){
	arma::mat Lchol;
	arma::mat U;
	arma::vec s;
	arma::mat V;

	arma::mat B=(RHO + RHO.t())/2.0;
	arma::svd(U,s,V,B);
	arma::mat H=V*arma::diagmat(s)*V.t();
	arma::mat RHO_hat=(B + H)/2.0;
	RHO_hat=(RHO_hat + RHO_hat.t())/2.0;

	bool ok=0;
	int k=0;
	while (ok==0){
		ok = arma::chol(Lchol, RHO_hat);
		k=k+1;
		if(ok==0){
			double mineig=arma::min(arma::eig_sym(RHO_hat));
			RHO_hat=RHO_hat + (-mineig*k*k + arma::eps(mineig))*arma::eye(arma::size(RHO));
		}

	}
	RHO=RHO_hat;

	return;
}

///////////////////////////////////
// A helper function for boost program opptions
///////////////////////////////////
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v){
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}

#endif
