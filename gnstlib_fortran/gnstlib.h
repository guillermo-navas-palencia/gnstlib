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

#include <complex>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif 

/* Exponential and logarithmic functions */
/******************************************************************************/
/* compute exp(x) */
double gnstlib_exp(const double x);
/* compute exp2(x) = 2^x */
double gnstlib_exp2(const double x);
/* compute exp10(x) = 10^x */
double gnstlib_exp10(const double x);
/* compute exp(x) - 1 */
double gnstlib_expm1(const double x);
/* compute (exp(x) - 1) / x */
double gnstlib_expm1x(const double x);
/* compute (exp(x) - 1 - x) / (0.5 * x * x) */
double gnstlib_expm1mx(const double x);

/* compute log(x) */
double gnstlib_log(const double x);
/* compute log2(x) */
double gnstlib_log2(const double x);
/* compute log10(x) */
double gnstlib_log10(const double x);
/* compute log(1 + x) */
double gnstlib_log1p(const double x);
/* compute log(1 + x) - x = log1p(x) - x */
double gnstlib_log1pmx(const double x); 
/* compute log(1 - exp(x)) */
double gnstlib_log1mexp(const double x, int& err_id);
/* compute log(1 + exp(x)) */
double gnstlib_log1pexp(const double x);

/* compute exp(z) */
std::complex<double> gnstlib_cexp(const std::complex<double> z);
/* compute log(z) */
std::complex<double> gnstlib_clog(const std::complex<double> z);

/* Trigonometric functions */
/******************************************************************************/
// compute cos(x)
double gnstlib_cos(const double x);
// compute sin(x)
double gnstlib_sin(const double x);
// compute tan(x)
double gnstlib_tan(const double x);
// compute sec(x)
double gnstlib_sec(const double x);
// compute csc(x)
double gnstlib_csc(const double x);
// compute cot(x)
double gnstlib_cot(const double x);

// compute acos(x)
double gnstlib_acos(const double x);
// compute asin(x)
double gnstlib_asin(const double x);
// compute atan(x)
double gnstlib_atan(const double x);
// còmpute atan2(y, x)
double gnstlib_atan2(const double y, const double x);
// compute asec(x)
double gnstlib_asec(const double x);
// compute acsc(x)
double gnstlib_acsc(const double x);
// compute acot(x)
double gnstlib_acot(const double x);

// compute cos(pi * x)
double gnstlib_cospi(const double x);
// compute sin(pi * x)
double gnstlib_sinpi(const double x);

// compute cos(z)
std::complex<double> gnstlib_ccos(const std::complex<double> z);
// compute sin(z)
std::complex<double> gnstlib_csin(const std::complex<double> z);
// compute tan(z)
std::complex<double> gnstlib_ctan(const std::complex<double> z);
// compute sec(z)
std::complex<double> gnstlib_csec(const std::complex<double> z);
// compute csc(z)
std::complex<double> gnstlib_ccsc(const std::complex<double> z);
// compute cot(z)
std::complex<double> gnstlib_ccot(const std::complex<double> z);

// compute acos(z)
std::complex<double> gnstlib_cacos(const std::complex<double> z);
// compute asin(z)
std::complex<double> gnstlib_casin(const std::complex<double> z);
// compute atan(z)
std::complex<double> gnstlib_catan(const std::complex<double> z);
// compute asec(z)
std::complex<double> gnstlib_casec(const std::complex<double> z);
// compute acsc(z)
std::complex<double> gnstlib_cacsc(const std::complex<double> z);
// compute acot(z)
std::complex<double> gnstlib_cacot(const std::complex<double> z);

// compute cos(pi * z)
std::complex<double> gnstlib_ccospi(const std::complex<double> z);
// compute sin(pi * z)
std::complex<double> gnstlib_csinpi(const std::complex<double> z);

/* Hyperbolic functions */
/******************************************************************************/
// compute cosh(x)
double gnstlib_cosh(const double x);
// compute sinh(x)
double gnstlib_sinh(const double x);
// compute tanh(x)
double gnstlib_tanh(const double x);
// compute sech(x)
double gnstlib_sech(const double x);
// compute csch(x)
double gnstlib_csch(const double x);
// compute coth(x)
double gnstlib_coth(const double x);

// compute acosh(x)
double gnstlib_acosh(const double x);
// compute asinh(x)
double gnstlib_asinh(const double x);
// compute atanh(x)
double gnstlib_atanh(const double x);
// compute asech(x)
double gnstlib_asech(const double x);
// compute acsch(x)
double gnstlib_acsch(const double x);
// compute acoth(x)
double gnstlib_acoth(const double x);

// compute cosh(z)
std::complex<double> gnstlib_ccosh(const std::complex<double> z);
// compute sinh(z)
std::complex<double> gnstlib_csinh(const std::complex<double> z);
// compute tanh(z)
std::complex<double> gnstlib_ctanh(const std::complex<double> z);
// compute sech(z)
std::complex<double> gnstlib_csech(const std::complex<double> z);
// compute csch(z)
std::complex<double> gnstlib_ccsch(const std::complex<double> z);
// compute coth(z)
std::complex<double> gnstlib_ccoth(const std::complex<double> z);

// compute acosh(z)
std::complex<double> gnstlib_cacosh(const std::complex<double> z);
// compute asinh(z)
std::complex<double> gnstlib_casinh(const std::complex<double> z);
// compute atanh(z)
std::complex<double> gnstlib_catanh(const std::complex<double> z);
// compute asech(z)
std::complex<double> gnstlib_casech(const std::complex<double> z);
// compute acsch(z)
std::complex<double> gnstlib_cacsch(const std::complex<double> z);
// compute acoth(z)
std::complex<double> gnstlib_cacoth(const std::complex<double> z);

/* Power functions */
/******************************************************************************/
// compute pow(x, y)
double gnstlib_pow(const double x, const double y);
// compute sqrt(x)
double gnstlib_sqrt(const double x);
// compute hypot(x, y)
double gnstlib_hypot(const double x, const double y);
// compute cbrt(x)
double gnstlib_cbrt(const double x);
// compute powm1(x, y) = x^y - 1
double gnstlib_powm1(const double x, const double y, int& err_id);

// compute power(cx, cy)
std::complex<double> gnstlib_cpow(const std::complex<double> x, 
                                  const std::complex<double> y);
// compute sqrt(z)
std::complex<double> gnstlib_csqrt(const std::complex<double> z);

/* Gamma functions */
/******************************************************************************/

/* main functions */

// compute gamma(x)
double gnstlib_gamma(const double x, int& err_id);
// compute gammaln(x) = Log(gamma(|x|))
double gnstlib_gammaln(const double x, int& err_id);
// compute gammasign(x) => gammaln(x)
int gnstlib_gammasign(const double x);
// compute factorial(n)
double gnstlib_factorial(const int n, int& err_id);
// compute quotient of two gammas
double gnstlib_qgamma(const double x, const double y, int& err_id);

// compute gamma(z)
std::complex<double> gnstlib_cgamma(const std::complex<double> z, 
                                    int& err_id);
// compute logGamma(z)
std::complex<double> gnstlib_loggamma(const std::complex<double> z, 
                                      int& err_id);
// compute rgamma(z) = 1 / gamma(z)
std::complex<double> gnstlib_rgamma(const std::complex<double> z);

/* auxiliary function */

// compute Stirling series (x > 0)
double gnstlib_stirling(const double x, int& err_id);
// compute scaled gamma function (x > 0)
double gnstlib_gammastar(const double x, int& err_id);
// compute auxiliary function g(x) for x in [-1,1]
double gnstlib_auxgam(const double x, int& err_id);

/* vectorized functions */

void gnstlib_gamma_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_gammaln_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_cgamma_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_loggamma_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_rgamma_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

/* Exponential, logarithmic and trigonometric integrals */
/******************************************************************************/

/* main functions */

// compute exponential integral ei(x)
double gnstlib_ei(const double x, int& err_id);
// compute exponential integral e1(x)
double gnstlib_e1(const double x, int& err_id);
// compute logarithmic integral li(x)
double gnstlib_li(const double x, int& err_id);
// compute cosine integral si(x)
double gnstlib_ci(const double x, int& err_id);
// compute sine integral si(x)
double gnstlib_si(const double x, int& err_id);
// compute hyperbolic cosine integral chi(x)
double gnstlib_chi(const double x, int& err_id);
// compute hyperbolic sine integral shi(x)
double gnstlib_shi(const double x, int& err_id);

// compute ei(z)
std::complex<double> gnstlib_cei(const std::complex<double> z, int& err_id);
// compute e1(z)
std::complex<double> gnstlib_ce1(const std::complex<double> z, int& err_id);
// compute li(z)
std::complex<double> gnstlib_cli(const std::complex<double> z, int& err_id);
// compute ci(z)
std::complex<double> gnstlib_cci(const std::complex<double> z, int& err_id);
// compute si(z)
std::complex<double> gnstlib_csi(const std::complex<double> z, int& err_id);
// compute chi(z)
std::complex<double> gnstlib_cchi(const std::complex<double> z, int& err_id);
// compute shi(z)
std::complex<double> gnstlib_cshi(const std::complex<double> z, int& err_id);

/* inverse functions */

// compute inverse of ei(x)
double gnstlib_invei(const double x, int& err_id);
// compute inverse of li(x)
double gnstlib_invli(const double x, int& err_id);

/* vectorized versions */

void gnstlib_ei_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_e1_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_li_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_ci_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_si_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_chi_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_shi_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_cei_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_ce1_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_cli_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_cci_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_csi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_cchi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_cshi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

/* Error functions, Dawson's and Fresnel Integrals */
/******************************************************************************/

/* main functions */

// compute erf(x)
double gnstlib_erf(const double x, int& err_id);
// compute erfc(x)
double gnstlib_erfc(const double x, int& err_id);
// compute erfcx(x)
double gnstlib_erfcx(const double x, int& err_id);
// compute erfi(x)
double gnstlib_erfi(const double x, int& err_id);
// compute dawson(x)
double gnstlib_dawson(const double x, int& err_id);
// compute fresnel C(x)
double gnstlib_fresnelc(const double x, int& err_id);
// compute fresnel S(x)
double gnstlib_fresnels(const double x, int& err_id);
// compute Voigt's profile
double gnstlib_voigt_profile(const double x, const double sigma, 
      const double gamma, int& err_id);

// compute erf(z)
std::complex<double> gnstlib_cerf(const std::complex<double> z, int& err_id);
// compute erfc(z)
std::complex<double> gnstlib_cerfc(const std::complex<double> z, int& err_id);
// compute erfcx(z)
std::complex<double> gnstlib_cerfcx(const std::complex<double> z, int& err_id);
// compute erfi(z)
std::complex<double> gnstlib_cerfi(const std::complex<double> z, int& err_id);
// compute dawson(z)
std::complex<double> gnstlib_cdawson(const std::complex<double> z, int& err_id);
// compute faddeeva(z)
std::complex<double> gnstlib_faddeeva(const std::complex<double> z, 
      int& err_id);
// compute fresnelc(z)
std::complex<double> gnstlib_cfresnelc(const std::complex<double> z, 
      int& err_id);
// compute fresnels(z)
std::complex<double> gnstlib_cfresnels(const std::complex<double> z, 
      int& err_id);

/* inverse functions */

// compute inverse erf(x)
double gnstlib_inverf(const double x, int& err_id);
// compute inverse erfc(x)
double gnstlib_inverfc(const double x, int& err_id);

/* vectorized versions */

void gnstlib_erf_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_erfc_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_erfcx_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_erfi_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_dawson_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_fresnelc_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_fresnels_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_inverf_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_inverfc_vec(const int n, const double *v, double *r, 
      unsigned short option = 0);

void gnstlib_cerf_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_cerfc_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_cerfcx_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_cerfi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_cdawson_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_faddeeva_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_cfresnelc_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

void gnstlib_cfresnels_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option = 0);

/* Incomplete gamma and generalized exponential integral */
/******************************************************************************/

/* main functions */

// compute regularized incomplete gamma function P(a, x)
double gnstlib_gammainc_p(const double a, const double x, int& err_id);
// compute regularized incomplete gamma function Q(a, x)
double gnstlib_gammainc_q(const double a, const double x, int& err_id);

// compute generalized exponential integral
double gnstlib_expint(const double v, const double x, int& err_id);
double gnstlib_expint_acc(const int n, const double e, const double x, 
  int& err_id);

#ifdef __cplusplus
}
#endif