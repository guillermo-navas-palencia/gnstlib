/* Copyright (C) 2017 A. Gil, G. Navas-Palencia, J. Segura and N. M. Temme
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

/**
 * gnstlib.hpp
 *  + mathematical constants 
 *  + elementary functions and vector operations
 *  + special function and vectorized versions
 */

#ifndef GNSTLIB_H
#define GNSTLIB_H

#include <complex>
#include <vector>

#include "gnstlib_basic.hpp"
#include "gnstlib_constants.hpp"

namespace GNSTLIB 
{

/* Gamma functions */
/******************************************************************************/

/* main functions */
  
// compute gamma(x)
double gamma(const double x, int& err_id);
// compute gammaln(x) = Log(gamma(|x|))
double gammaln(const double x, int& err_id);
// compute gammasign(x) => gammaln(x)
int gammasign(const double x);
// compute factorial(n)
double factorial(const int n, int& err_id);
// compute quotient of two gammas
double qgamma(const double x, const double y, int& err_id);

// compute gamma(z)
std::complex<double> gamma(const std::complex<double> z, int& err_id);
// compute logGamma(z)
std::complex<double> loggamma(const std::complex<double> z, int& err_id);
// compute rgamma(z) = 1 / gamma(z)
std::complex<double> rgamma(const std::complex<double> z);

/* auxiliary function */

// compute Stirling series (x > 0)
double stirling(const double x, int& err_id);
// compute scaled gamma function (x > 0)
double gammastar(const double x, int& err_id);
// compute auxiliary function g(x) for x in [-1,1]
double auxgam(const double x, int& err_id);

/* vectorized functions */

void gamma_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void gamma_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gammaln_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void gammaln_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gamma_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void gamma_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void loggamma_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void loggamma_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void rgamma_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void rgamma_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

/* Exponential, logarithmic and trigonometric integrals */
/******************************************************************************/

/* main functions */

// compute exponential integral ei(x)
double ei(const double x, int& err_id);
// compute exponential integral e1(x)
double e1(const double x, int& err_id);
// compute logarithmic integral li(x)
double li(const double x, int& err_id);
// compute cosine integral si(x)
double ci(const double x, int& err_id);
// compute sine integral si(x)
double si(const double x, int& err_id);
// compute hyperbolic cosine integral chi(x)
double chi(const double x, int& err_id);
// compute hyperbolic sine integral shi(x)
double shi(const double x, int& err_id);

// compute ei(z)
std::complex<double> ei(const std::complex<double> z, int& err_id);
// compute e1(z)
std::complex<double> e1(const std::complex<double> z, int& err_id);
// compute li(z)
std::complex<double> li(const std::complex<double> z, int& err_id);
// compute ci(z)
std::complex<double> ci(const std::complex<double> z, int& err_id);
// compute si(z)
std::complex<double> si(const std::complex<double> z, int& err_id);
// compute chi(z)
std::complex<double> chi(const std::complex<double> z, int& err_id);
// compute shi(z)
std::complex<double> shi(const std::complex<double> z, int& err_id);

/* extra functions */

void sici(const std::complex<double> z, std::complex<double>& si, 
      std::complex<double>& ci, int& err_id);

/* inverse functions */

// compute inverse of ei(x)
double invei(const double x, int& err_id);
// compute inverse of li(x)
double invli(const double x, int& err_id);

/* vectorized versions */

void ei_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void ei_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void e1_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void e1_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void li_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void li_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void ci_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void ci_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void si_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void si_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void chi_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void chi_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void shi_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void shi_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void ei_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void ei_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void e1_vec(std::vector<std::complex<double>>& v,
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void e1_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void li_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void li_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void ci_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void ci_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void si_vec(std::vector<std::complex<double>>& v,
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void si_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void chi_vec(std::vector<std::complex<double>>& v,
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void chi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void shi_vec(std::vector<std::complex<double>>& v,
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void shi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

/* Error functions, Dawson's and Fresnel Integrals */
/******************************************************************************/

/* main functions */

// compute erf(x)
double erf(const double x, int& err_id);
// compute erfc(x)
double erfc(const double x, int& err_id);
// compute erfcx(x)
double erfcx(const double x, int& err_id);
// compute erfi(x)
double erfi(const double x, int& err_id);
// compute dawson(x)
double dawson(const double x, int& err_id);
// compute fresnel C(x)
double fresnelc(const double x, int& err_id);
// compute fresnel S(x)
double fresnels(const double x, int& err_id);
// compute Voigt's profile
double voigt_profile(const double x, const double sigma, const double gamma, 
      int& err_id);

// compute erf(z)
std::complex<double> erf(const std::complex<double> z, int& err_id);
// compute erfc(z)
std::complex<double> erfc(const std::complex<double> z, int& err_id);
// compute erfcx(z)
std::complex<double> erfcx(const std::complex<double> z, int& err_id);
// compute erfi(z)
std::complex<double> erfi(const std::complex<double> z, int& err_id);
// compute dawson(z)
std::complex<double> dawson(const std::complex<double> z, int& err_id);
// compute faddeeva(z)
std::complex<double> faddeeva(const std::complex<double> z, int& err_id);
// compute fresnelc(z)
std::complex<double> fresnelc(const std::complex<double> z, int& err_id);
// compute fresnels(z)
std::complex<double> fresnels(const std::complex<double> z, int& err_id);

/* inverse functions */

// compute inverse erf(x)
double inverf(const double x, int& err_id);
// compute inverse erfc(x)
double inverfc(const double x, int& err_id);

/* vectorized versions */

void erf_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void erf_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void erfc_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void erfc_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void erfcx_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void erfcx_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void erfi_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void erfi_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void dawson_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void dawson_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void fresnelc_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void fresnelc_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void fresnels_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void fresnels_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void inverf_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void inverf_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void inverfc_vec(std::vector<double>& v, std::vector<double>& r, 
      unsigned short option = 0);
void inverfc_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void erf_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void erf_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void erfc_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void erfc_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void erfcx_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void erfcx_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void erfi_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void erfi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void dawson_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void dawson_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void faddeeva_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void faddeeva_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void fresnelc_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void fresnelc_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void fresnels_vec(std::vector<std::complex<double>>& v, 
      std::vector<std::complex<double>>& r, unsigned short option = 0);
void fresnels_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

/* Incomplete gamma and generalized exponential integral */
/******************************************************************************/

/* main functions */

// compute regularized incomplete gamma function P(a, x)
double gammainc_p(const double a, const double x, int& err_id);
// compute regularized incomplete gamma function Q(a, x)
double gammainc_q(const double a, const double x, int& err_id);

// compute generalized exponential integral
double expint(const double v, const double x, int& err_id);
double expint(const int n, const double e, const double x, int& err_id);

} // namespace GNSTLIB

#endif /* GNSTLIB_H */
