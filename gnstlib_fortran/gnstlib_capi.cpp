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
#include <iterator>
#include <vector>

#include "gnstlib.h"
#include "gnstlib.hpp"

// Exponential and logarithmic functions

double gnstlib_exp(const double x)
{
  return GNSTLIB::exp(x);
}

double gnstlib_exp2(const double x)
{
  return GNSTLIB::exp2(x);
}

double gnstlib_exp10(const double x)
{
  return GNSTLIB::exp10(x);
}

double gnstlib_expm1(const double x)
{
  return GNSTLIB::expm1(x);
}

double gnstlib_expm1x(const double x)
{
  return GNSTLIB::expm1x(x);
}

double gnstlib_expm1mx(const double x)
{
  return GNSTLIB::expm1mx(x);
}

double gnstlib_log(const double x)
{
  return GNSTLIB::log(x);
}

double gnstlib_log2(const double x)
{
  return GNSTLIB::log2(x);
}

double gnstlib_log10(const double x)
{
  return GNSTLIB::log10(x);
}

double gnstlib_log1p(const double x)
{
  return GNSTLIB::log1p(x);
}

double gnstlib_log1pmx(const double x)
{
  return GNSTLIB::log1pmx(x);
}

double gnstlib_log1mexp(const double x, int& err_id)
{
  return GNSTLIB::log1mexp(x, err_id);
}

double gnstlib_log1pexp(const double x)
{
  return GNSTLIB::log1pexp(x);
}

std::complex<double> gnstlib_cexp(const std::complex<double> z)
{
  return GNSTLIB::exp(z);
}

std::complex<double> gnstlib_clog(const std::complex<double> z)
{
  return GNSTLIB::log(z);
}

// Trigonometric functions

double gnstlib_cos(const double x)
{
  return GNSTLIB::cos(x);
}

double gnstlib_sin(const double x)
{
  return GNSTLIB::sin(x);
}

double gnstlib_tan(const double x)
{
  return GNSTLIB::tan(x);
}

double gnstlib_sec(const double x)
{
  return GNSTLIB::sec(x);
}

double gnstlib_csc(const double x)
{
  return GNSTLIB::csc(x);
}

double gnstlib_cot(const double x)
{
  return GNSTLIB::cot(x);
}

double gnstlib_acos(const double x)
{
  return GNSTLIB::acos(x);
}

double gnstlib_asin(const double x)
{
  return GNSTLIB::asin(x);
}

double gnstlib_atan(const double x)
{
  return GNSTLIB::atan(x);
}

double gnstlib_atan2(const double y, const double x)
{
  return GNSTLIB::atan2(y, x);
}

double gnstlib_asec(const double x)
{
  return GNSTLIB::asec(x);
}

double gnstlib_acsc(const double x)
{
  return GNSTLIB::acsc(x);
}

double gnstlib_acot(const double x)
{
  return GNSTLIB::acot(x);
}

double gnstlib_cospi(const double x)
{
  return GNSTLIB::cospi(x);
}

double gnstlib_sinpi(const double x)
{
  return GNSTLIB::sinpi(x);
}

std::complex<double> gnstlib_ccos(const std::complex<double> z)
{
  return GNSTLIB::cos(z);
}

std::complex<double> gnstlib_csin(const std::complex<double> z)
{
  return GNSTLIB::sin(z);
}

std::complex<double> gnstlib_ctan(const std::complex<double> z)
{
  return GNSTLIB::tan(z);
}

std::complex<double> gnstlib_csec(const std::complex<double> z)
{
  return GNSTLIB::sec(z);
}

std::complex<double> gnstlib_ccsc(const std::complex<double> z)
{
  return GNSTLIB::csc(z);
}

std::complex<double> gnstlib_ccot(const std::complex<double> z)
{
  return GNSTLIB::cot(z);
}

std::complex<double> gnstlib_cacos(const std::complex<double> z)
{
  return GNSTLIB::acos(z);
}

std::complex<double> gnstlib_casin(const std::complex<double> z)
{
  return GNSTLIB::asin(z);
}

std::complex<double> gnstlib_catan(const std::complex<double> z)
{
  return GNSTLIB::atan(z);
}

std::complex<double> gnstlib_casec(const std::complex<double> z)
{
  return GNSTLIB::asec(z);
}

std::complex<double> gnstlib_cacsc(const std::complex<double> z)
{
  return GNSTLIB::acsc(z);
}

std::complex<double> gnstlib_cacot(const std::complex<double> z)
{
  return GNSTLIB::acot(z);
}

std::complex<double> gnstlib_ccospi(const std::complex<double> z)
{
  return GNSTLIB::cospi(z);
}

std::complex<double> gnstlib_csinpi(const std::complex<double> z)
{
  return GNSTLIB::sinpi(z);
}

// Hyperbolic functions

double gnstlib_cosh(const double x)
{
  return GNSTLIB::cosh(x);
}

double gnstlib_sinh(const double x)
{
  return GNSTLIB::sinh(x);
}

double gnstlib_tanh(const double x)
{
  return GNSTLIB::tanh(x);
}

double gnstlib_sech(const double x)
{
  return GNSTLIB::sech(x);
}

double gnstlib_csch(const double x)
{
  return GNSTLIB::csch(x);
}

double gnstlib_coth(const double x)
{
  return GNSTLIB::coth(x);
}

double gnstlib_acosh(const double x)
{
  return GNSTLIB::acosh(x);
}

double gnstlib_asinh(const double x)
{
  return GNSTLIB::asinh(x);
}

double gnstlib_atanh(const double x)
{
  return GNSTLIB::atanh(x);
}

double gnstlib_asech(const double x)
{
  return GNSTLIB::asech(x);
}

double gnstlib_acsch(const double x)
{
  return GNSTLIB::acsch(x);
}

double gnstlib_acoth(const double x)
{
  return GNSTLIB::acoth(x);
}

std::complex<double> gnstlib_ccosh(const std::complex<double> z)
{
  return GNSTLIB::cosh(z);
}

std::complex<double> gnstlib_csinh(const std::complex<double> z)
{
  return GNSTLIB::sinh(z);
}

std::complex<double> gnstlib_ctanh(const std::complex<double> z)
{
  return GNSTLIB::tanh(z);
}

std::complex<double> gnstlib_csech(const std::complex<double> z)
{
  return GNSTLIB::sech(z);
}

std::complex<double> gnstlib_ccsch(const std::complex<double> z)
{
  return GNSTLIB::csch(z);
}

std::complex<double> gnstlib_ccoth(const std::complex<double> z)
{
  return GNSTLIB::coth(z);
}

std::complex<double> gnstlib_cacosh(const std::complex<double> z)
{
  return GNSTLIB::acosh(z);
}

std::complex<double> gnstlib_casinh(const std::complex<double> z)
{
  return GNSTLIB::asinh(z);
}

std::complex<double> gnstlib_catanh(const std::complex<double> z)
{
  return GNSTLIB::atanh(z);
}

std::complex<double> gnstlib_casech(const std::complex<double> z)
{
  return GNSTLIB::asech(z);
}

std::complex<double> gnstlib_cacsch(const std::complex<double> z)
{
  return GNSTLIB::acsch(z);
}

std::complex<double> gnstlib_cacoth(const std::complex<double> z)
{
  return GNSTLIB::acoth(z);
}

// Power functions

double gnstlib_pow(const double x, const double y)
{
  return GNSTLIB::pow(x, y);
}

double gnstlib_sqrt(const double x)
{
  return GNSTLIB::sqrt(x);
}

double gnstlib_hypot(const double x, const double y)
{
  return GNSTLIB::hypot(x, y);
}

double gnstlib_cbrt(const double x)
{
  return GNSTLIB::cbrt(x);
}

double gnstlib_powm1(const double x, const double y, int& err_id)
{
  return GNSTLIB::powm1(x, y, err_id);
}

std::complex<double> gnstlib_cpow(const std::complex<double> x, 
                                  const std::complex<double> y)
{
  return GNSTLIB::pow(x, y);
}

std::complex<double> gnstlib_csqrt(const std::complex<double> z)
{
  return GNSTLIB::sqrt(z);
}

// Gamma functions

double gnstlib_gamma(const double x, int& err_id)
{
  return GNSTLIB::gamma(x, err_id);
}

double gnstlib_gammaln(const double x, int& err_id)
{
  return GNSTLIB::gammaln(x, err_id);
}

int gnstlib_gammasign(const double x)
{
  return GNSTLIB::gammasign(x);
}

double gnstlib_factorial(const int n, int& err_id)
{
  return GNSTLIB::factorial(n, err_id);
}

double gnstlib_qgamma(const double x, const double y, int& err_id)
{
  return GNSTLIB::qgamma(x, y, err_id);
}

std::complex<double> gnstlib_cgamma(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::gamma(z, err_id);
}

std::complex<double> gnstlib_loggamma(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::loggamma(z, err_id);
}

std::complex<double> gnstlib_rgamma(const std::complex<double> z)
{
  return GNSTLIB::rgamma(z);
}

double gnstlib_stirling(const double x, int& err_id)
{
  return GNSTLIB::stirling(x, err_id);
}

double gnstlib_gammastar(const double x, int& err_id)
{
  return GNSTLIB::gammastar(x, err_id);
}

double gnstlib_auxgam(const double x, int& err_id)
{
  return GNSTLIB::auxgam(x, err_id);
}

void gnstlib_gamma_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::gamma_vec(n, v, r, option);
}

void gnstlib_gammaln_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::gammaln_vec(n, v, r, option);
}

void gnstlib_cgamma_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::gamma_vec(n, v, r, option);
}

void gnstlib_loggamma_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::loggamma_vec(n, v, r, option);
}

void gnstlib_rgamma_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::rgamma_vec(n, v, r, option);
}

// Exponential, logarithmic and trigonometric integrals

double gnstlib_ei(const double x, int& err_id)
{
  return GNSTLIB::ei(x, err_id);
}

double gnstlib_e1(const double x, int& err_id)
{
  return GNSTLIB::e1(x, err_id);
}

double gnstlib_li(const double x, int& err_id)
{
  return GNSTLIB::li(x, err_id);
}

double gnstlib_ci(const double x, int& err_id)
{
  return GNSTLIB::ci(x, err_id);
}

double gnstlib_si(const double x, int& err_id)
{
  return GNSTLIB::si(x, err_id);
}

double gnstlib_chi(const double x, int& err_id)
{
  return GNSTLIB::chi(x, err_id);
}

double gnstlib_shi(const double x, int& err_id)
{
  return GNSTLIB::shi(x, err_id);
}

std::complex<double> gnstlib_cei(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::ei(z, err_id);
}

std::complex<double> gnstlib_ce1(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::e1(z, err_id);
}

std::complex<double> gnstlib_cli(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::li(z, err_id);
}

std::complex<double> gnstlib_cci(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::ci(z, err_id);
}

std::complex<double> gnstlib_csi(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::si(z, err_id);
}

std::complex<double> gnstlib_cchi(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::chi(z, err_id);
}

std::complex<double> gnstlib_cshi(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::shi(z, err_id);
}

double gnstlib_invei(const double x, int& err_id)
{
  return GNSTLIB::invei(x, err_id);
}

double gnstlib_invli(const double x, int& err_id)
{
  return GNSTLIB::invei(x, err_id);
}

void gnstlib_ei_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::ei_vec(n, v, r, option);
}

void gnstlib_e1_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::e1_vec(n, v, r, option);
}

void gnstlib_li_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::li_vec(n, v, r, option);
}

void gnstlib_ci_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::ci_vec(n, v, r, option);
}

void gnstlib_si_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::si_vec(n, v, r, option);
}

void gnstlib_chi_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::chi_vec(n, v, r, option);
}

void gnstlib_shi_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::shi_vec(n, v, r, option);
}

void gnstlib_cei_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::ei_vec(n, v, r, option);
}

void gnstlib_ce1_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::e1_vec(n, v, r, option);
}

void gnstlib_cli_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::li_vec(n, v, r, option);
}

void gnstlib_cci_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::ci_vec(n, v, r, option);
}

void gnstlib_csi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::si_vec(n, v, r, option);
}

void gnstlib_cchi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::chi_vec(n, v, r, option);
}

void gnstlib_cshi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::shi_vec(n, v, r, option);
}

// Error functions, Dawson's and Fresnel Integrals

double gnstlib_erf(const double x, int& err_id)
{
  return GNSTLIB::erf(x, err_id);
}

double gnstlib_erfc(const double x, int& err_id)
{
  return GNSTLIB::erfc(x, err_id);
}

double gnstlib_erfcx(const double x, int& err_id)
{
  return GNSTLIB::erfcx(x, err_id);
}

double gnstlib_erfi(const double x, int& err_id)
{
  return GNSTLIB::erfi(x, err_id);
}

double gnstlib_dawson(const double x, int& err_id)
{
  return GNSTLIB::dawson(x, err_id);
}

double gnstlib_fresnelc(const double x, int& err_id)
{
  return GNSTLIB::fresnelc(x, err_id);
}

double gnstlib_fresnels(const double x, int& err_id)
{
  return GNSTLIB::fresnels(x, err_id);
}

double gnstlib_voigt_profile(const double x, const double sigma, 
  const double gamma, int& err_id)
{
  return GNSTLIB::voigt_profile(x, sigma, gamma, err_id); 
}

std::complex<double> gnstlib_cerf(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::erf(z, err_id);
}

std::complex<double> gnstlib_cerfc(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::erfc(z, err_id);
}

std::complex<double> gnstlib_cerfcx(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::erfcx(z, err_id);
}

std::complex<double> gnstlib_cerfi(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::erfi(z, err_id);
}

std::complex<double> gnstlib_cdawson(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::dawson(z, err_id);
}

std::complex<double> gnstlib_faddeeva(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::faddeeva(z, err_id);
}

std::complex<double> gnstlib_cfresnelc(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::fresnelc(z, err_id);
}

std::complex<double> gnstlib_cfresnels(const std::complex<double> z, int& err_id)
{
  return GNSTLIB::fresnels(z, err_id);
}

double gnstlib_inverf(const double x, int& err_id)
{
  return GNSTLIB::inverf(x, err_id);
}

double gnstlib_inverfc(const double x, int& err_id)
{
  return GNSTLIB::inverfc(x, err_id);
}

void gnstlib_erf_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::erf_vec(n, v, r, option);
}

void gnstlib_erfc_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::erfc_vec(n, v, r, option);
}

void gnstlib_erfcx_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::erfcx_vec(n, v, r, option);
}

void gnstlib_erfi_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::erfi_vec(n, v, r, option);
}

void gnstlib_dawson_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::dawson_vec(n, v, r, option);
}

void gnstlib_fresnelc_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::fresnelc_vec(n, v, r, option);
}

void gnstlib_fresnels_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::fresnels_vec(n, v, r, option);
}

void gnstlib_inverf_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::inverf_vec(n, v, r, option);
}

void gnstlib_inverfc_vec(const int n, const double *v, double *r, 
      unsigned short option)
{
  GNSTLIB::inverfc_vec(n, v, r, option);
}

void gnstlib_cerf_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::erf_vec(n, v, r, option);
}

void gnstlib_cerfc_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::erfc_vec(n, v, r, option);
}

void gnstlib_cerfcx_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::erfcx_vec(n, v, r, option);
}

void gnstlib_cerfi_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::erfi_vec(n, v, r, option);
}

void gnstlib_cdawson_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::dawson_vec(n, v, r, option);
}

void gnstlib_faddeeva_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::faddeeva_vec(n, v, r, option);
}

void gnstlib_cfresnelc_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::fresnelc_vec(n, v, r, option);
}

void gnstlib_cfresnels_vec(const int n, const std::complex<double> *v, 
      std::complex<double> *r, unsigned short option)
{
  GNSTLIB::fresnels_vec(n, v, r, option);
}

// Incomplete gamma and generalized exponential integral

double gnstlib_gammainc_p(const double a, const double x, int& err_id)
{
  return GNSTLIB::gammainc_p(a, x, err_id);
}

double gnstlib_gammainc_q(const double a, const double x, int& err_id)
{
  return GNSTLIB::gammainc_q(a, x, err_id);
}

double gnstlib_expint(const double v, const double x, int& err_id)
{
  return GNSTLIB::expint(v, x, err_id);
}

double gnstlib_expint_acc(const int n, const double e, const double x, 
  int& err_id)
{
  return GNSTLIB::expint(n, e, x, err_id);
}  