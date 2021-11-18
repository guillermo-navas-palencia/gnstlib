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

#include <cmath>
#include <functional>
#include <vector>
#include <omp.h>

#include "gnstlib.hpp"

typedef std::complex<double> cmplx;

// vectorize univariate special function - serial mode
template<typename F, typename T>
void vectorize_seq(std::function<F(T)> f, std::vector<T>& v, std::vector<F>& r) 
{
  for (auto i = begin(v); i != end(v); ++i)
    r.push_back(f(*i));
}

template<typename F, typename T>
void vectorize_seq(std::function<F(T)> f, const int n, const T *v, F *r)
{
  for (int i = 0; i < n; i++)
    r[i] = f(v[i]);
}

// vectorize univariate special function - parallel mode
template<typename F, typename T>
void vectorize_omp(std::function<F(T)> f, std::vector<T>& v, std::vector<F>& r) 
{
  int thread;
  size_t i, sz;
  sz = v.size();

  // create private vector to fill unordered values of f(i) and return
  // ordered array r in parallel

  #pragma omp parallel
  {
    std::vector<F> r_priv;
    #pragma omp for nowait schedule(static)
    for (i = 0; i < sz; i++) {
      r_priv.push_back(f(v[i]));
    }

    #pragma omp for schedule(static) ordered
    for (thread = 0; thread < omp_get_num_threads(); thread++) {
      #pragma omp ordered
      r.insert(r.end(), r_priv.begin(), r_priv.end());
    }
  }
}

// vectorize univariate special function - parallel mode
template<typename F, typename T>
void vectorize_omp(std::function<F(T)> f, const int n, const T *v, F* r) 
{
  int i;

  // fixed size, not reallocation (thread-safe)
  #pragma omp parallel for
  for (i = 0; i < n; i++)
    r[i] = f(v[i]);
}

// vectorize univariate special functions
template<typename F, typename T>
void vectorize(std::function<F(T)> f, std::vector<T>& v, std::vector<F>& r, 
  const unsigned short pmode) 
{
  switch(pmode) {
    case 0: vectorize_seq(f, v, r); break;
    case 1: vectorize_omp(f, v, r); break;
    default: vectorize_seq(f, v, r); break;
  }
}

template<typename F, typename T>
void vectorize(std::function<F(T)> f, const int n, const T *v, F *r, 
  const unsigned short pmode) 
{
  switch(pmode) {
    case 0: vectorize_seq(f, n, v, r); break;
    case 1: vectorize_omp(f, n, v, r); break;
    default: vectorize_seq(f, n, v, r); break;
  }
}

// Special functions:
//  general comment: err_id is dummy to enhance performance.

/* Gamma functions */
/******************************************************************************/
void GNSTLIB::gamma_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_gamma = 
    [](double x) {int err_id; return GNSTLIB::gamma(x, err_id);};

  vectorize(f_gamma, v, r, option);
}

void GNSTLIB::gamma_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_gamma = 
    [](double x) {int err_id; return GNSTLIB::gamma(x, err_id);};

  vectorize(f_gamma, n, v, r, option);
}

void GNSTLIB::gammaln_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_gammaln = 
    [](double x) {int err_id; return GNSTLIB::gammaln(x, err_id);};

  vectorize(f_gammaln, v, r, option);
}

void GNSTLIB::gammaln_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_gammaln = 
    [](double x) {int err_id; return GNSTLIB::gamma(x, err_id);};

  vectorize(f_gammaln, n, v, r, option);
}

void GNSTLIB::gamma_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option) 
{
  std::function<cmplx(cmplx)> f_gamma = 
    [](cmplx z) {int err_id; return GNSTLIB::gamma(z, err_id);};

  vectorize(f_gamma, v, r, option);
}

void GNSTLIB::gamma_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_gamma = 
    [](cmplx z) {int err_id; return GNSTLIB::gamma(z, err_id);};

  vectorize(f_gamma, n, v, r, option);
}

void GNSTLIB::loggamma_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option) 
{
  std::function<cmplx(cmplx)> f_loggamma = 
    [](cmplx z) {int err_id; return GNSTLIB::loggamma(z, err_id);};

  vectorize(f_loggamma, v, r, option);
}

void GNSTLIB::loggamma_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_loggamma = 
    [](cmplx z) {int err_id; return GNSTLIB::loggamma(z, err_id);};

  vectorize(f_loggamma, n, v, r, option);
}

void GNSTLIB::rgamma_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option) 
{
  std::function<cmplx(cmplx)> f_rgamma = 
    [](cmplx z) {return GNSTLIB::rgamma(z);};

  vectorize(f_rgamma, v, r, option);
}

void GNSTLIB::rgamma_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_rgamma = 
    [](cmplx z) {return GNSTLIB::rgamma(z);};

  vectorize(f_rgamma, n, v, r, option);
}

/* Exponential, logarithmic and trigonometric integrals */
/******************************************************************************/
void GNSTLIB::ei_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_ei = 
    [](double x) {int err_id; return GNSTLIB::ei(x, err_id);};

  vectorize(f_ei, v, r, option);
}

void GNSTLIB::ei_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_ei = 
    [](double x) {int err_id; return GNSTLIB::ei(x, err_id);};

  vectorize(f_ei, n, v, r, option);
}

void GNSTLIB::ei_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_ei = 
    [](cmplx z) {int err_id; return GNSTLIB::ei(z, err_id);};

  vectorize(f_ei, v, r, option);
}

void GNSTLIB::ei_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_ei = 
    [](cmplx z) {int err_id; return GNSTLIB::ei(z, err_id);};

  vectorize(f_ei, n, v, r, option);
}

void GNSTLIB::e1_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_e1 = 
    [](double x) {int err_id; return GNSTLIB::e1(x, err_id);};

  vectorize(f_e1, v, r, option);
}

void GNSTLIB::e1_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_e1 = 
    [](double x) {int err_id; return GNSTLIB::e1(x, err_id);};

  vectorize(f_e1, n, v, r, option);
}

void GNSTLIB::e1_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_e1 = 
    [](cmplx z) {int err_id; return GNSTLIB::e1(z, err_id);};

  vectorize(f_e1, v, r, option);
}

void GNSTLIB::e1_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_e1 = 
    [](cmplx z) {int err_id; return GNSTLIB::e1(z, err_id);};

  vectorize(f_e1, n, v, r, option);
}

// void GNSTLIB::e1_vec_mkl(const int n, const double *v, double *r)
// {
//  vdExpInt1(n, v, r);
// }

void GNSTLIB::li_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_li = 
    [](double x) {int err_id; return GNSTLIB::li(x, err_id);};

  vectorize(f_li, v, r, option);
}

void GNSTLIB::li_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_li = 
    [](double x) {int err_id; return GNSTLIB::li(x, err_id);};

  vectorize(f_li, n, v, r, option);
}

void GNSTLIB::li_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_li =
    [](cmplx z) {int err_id; return GNSTLIB::li(z, err_id);};

  vectorize(f_li, v, r, option);
}

void GNSTLIB::li_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_li = 
    [](cmplx z) {int err_id; return GNSTLIB::li(z, err_id);};

  vectorize(f_li, n, v, r, option);
}

void GNSTLIB::ci_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_ci = 
    [](double x) {int err_id; return GNSTLIB::ci(x, err_id);};

  vectorize(f_ci, v, r, option);
}

void GNSTLIB::ci_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_ci = 
    [](double x) {int err_id; return GNSTLIB::ci(x, err_id);};

  vectorize(f_ci, n, v, r, option);
}

void GNSTLIB::ci_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_ci =
    [](cmplx z) {int err_id; return GNSTLIB::ci(z, err_id);};

  vectorize(f_ci, v, r, option);
}

void GNSTLIB::ci_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_ci = 
    [](cmplx z) {int err_id; return GNSTLIB::ci(z, err_id);};

  vectorize(f_ci, n, v, r, option);
}

void GNSTLIB::si_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_si = 
    [](double x) {int err_id; return GNSTLIB::si(x, err_id);};

  vectorize(f_si, v, r, option);
}

void GNSTLIB::si_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_si = 
    [](double x) {int err_id; return GNSTLIB::si(x, err_id);};

  vectorize(f_si, n, v, r, option);
}

void GNSTLIB::si_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_si =
    [](cmplx z) {int err_id; return GNSTLIB::si(z, err_id);};

  vectorize(f_si, v, r, option);
}

void GNSTLIB::si_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_si = 
    [](cmplx z) {int err_id; return GNSTLIB::si(z, err_id);};

  vectorize(f_si, n, v, r, option);
}

void GNSTLIB::chi_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_chi = 
    [](double x) {int err_id; return GNSTLIB::chi(x, err_id);};

  vectorize(f_chi, v, r, option);
}

void GNSTLIB::chi_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_chi = 
    [](double x) {int err_id; return GNSTLIB::chi(x, err_id);};

  vectorize(f_chi, n, v, r, option);
}

void GNSTLIB::chi_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_chi =
    [](cmplx z) {int err_id; return GNSTLIB::chi(z, err_id);};

  vectorize(f_chi, v, r, option);
}

void GNSTLIB::chi_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_chi = 
    [](cmplx z) {int err_id; return GNSTLIB::chi(z, err_id);};

  vectorize(f_chi, n, v, r, option);
}

void GNSTLIB::shi_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_shi = 
    [](double x) {int err_id; return GNSTLIB::shi(x, err_id);};

  vectorize(f_shi, v, r, option);
}

void GNSTLIB::shi_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_shi = 
    [](double x) {int err_id; return GNSTLIB::shi(x, err_id);};

  vectorize(f_shi, n, v, r, option);
}

void GNSTLIB::shi_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_shi =
    [](cmplx z) {int err_id; return GNSTLIB::shi(z, err_id);};

  vectorize(f_shi, v, r, option);
}

void GNSTLIB::shi_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_shi = 
    [](cmplx z) {int err_id; return GNSTLIB::shi(z, err_id);};

  vectorize(f_shi, n, v, r, option);
}

/* Error functions, Dawson's and Fresnel Integrals */
/******************************************************************************/

void GNSTLIB::erf_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_erf = 
    [](double x) {int err_id; return GNSTLIB::erf(x, err_id);};

  vectorize(f_erf, v, r, option);
}

void GNSTLIB::erf_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_erf = 
    [](double x) {int err_id; return GNSTLIB::erf(x, err_id);};

  vectorize(f_erf, n, v, r, option);
}

void GNSTLIB::erf_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_erf =
    [](cmplx z) {int err_id; return GNSTLIB::erf(z, err_id);};

  vectorize(f_erf, v, r, option);
}

void GNSTLIB::erf_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_erf = 
    [](cmplx z) {int err_id; return GNSTLIB::erf(z, err_id);};

  vectorize(f_erf, n, v, r, option);
}

void GNSTLIB::erfc_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_erfc = 
    [](double x) {int err_id; return GNSTLIB::erfc(x, err_id);};

  vectorize(f_erfc, v, r, option);
}

void GNSTLIB::erfc_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_erfc = 
    [](double x) {int err_id; return GNSTLIB::erfc(x, err_id);};

  vectorize(f_erfc, n, v, r, option);
}

void GNSTLIB::erfc_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_erfc =
    [](cmplx z) {int err_id; return GNSTLIB::erfc(z, err_id);};

  vectorize(f_erfc, v, r, option);
}

void GNSTLIB::erfc_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_erfc = 
    [](cmplx z) {int err_id; return GNSTLIB::erfc(z, err_id);};

  vectorize(f_erfc, n, v, r, option);
}

void GNSTLIB::erfcx_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_erfcx = 
    [](double x) {int err_id; return GNSTLIB::erfcx(x, err_id);};

  vectorize(f_erfcx, v, r, option);
}

void GNSTLIB::erfcx_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_erfcx = 
    [](double x) {int err_id; return GNSTLIB::erfcx(x, err_id);};

  vectorize(f_erfcx, n, v, r, option);
}

void GNSTLIB::erfcx_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_erfcx =
    [](cmplx z) {int err_id; return GNSTLIB::erfcx(z, err_id);};

  vectorize(f_erfcx, v, r, option);
}

void GNSTLIB::erfcx_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_erfcx = 
    [](cmplx z) {int err_id; return GNSTLIB::erfcx(z, err_id);};

  vectorize(f_erfcx, n, v, r, option);
}

void GNSTLIB::erfi_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_erfi = 
    [](double x) {int err_id; return GNSTLIB::erfi(x, err_id);};

  vectorize(f_erfi, v, r, option);
}

void GNSTLIB::erfi_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_erfi = 
    [](double x) {int err_id; return GNSTLIB::erfi(x, err_id);};

  vectorize(f_erfi, n, v, r, option);
}

void GNSTLIB::erfi_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_erfi =
    [](cmplx z) {int err_id; return GNSTLIB::erfi(z, err_id);};

  vectorize(f_erfi, v, r, option);
}

void GNSTLIB::erfi_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_erfi = 
    [](cmplx z) {int err_id; return GNSTLIB::erfi(z, err_id);};

  vectorize(f_erfi, n, v, r, option);
}

void GNSTLIB::dawson_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_dawson = 
    [](double x) {int err_id; return GNSTLIB::dawson(x, err_id);};

  vectorize(f_dawson, v, r, option);
}

void GNSTLIB::dawson_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_dawson = 
    [](double x) {int err_id; return GNSTLIB::dawson(x, err_id);};

  vectorize(f_dawson, n, v, r, option);
}

void GNSTLIB::dawson_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_dawson =
    [](cmplx z) {int err_id; return GNSTLIB::dawson(z, err_id);};

  vectorize(f_dawson, v, r, option);
}

void GNSTLIB::dawson_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_dawson = 
    [](cmplx z) {int err_id; return GNSTLIB::dawson(z, err_id);};

  vectorize(f_dawson, n, v, r, option);
}

void GNSTLIB::faddeeva_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_w =
    [](cmplx z) {int err_id; return GNSTLIB::faddeeva(z, err_id);};

  vectorize(f_w, v, r, option);
}

void GNSTLIB::faddeeva_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_w = 
    [](cmplx z) {int err_id; return GNSTLIB::faddeeva(z, err_id);};

  vectorize(f_w, n, v, r, option);
}

void GNSTLIB::fresnelc_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_fc = 
    [](double x) {int err_id; return GNSTLIB::fresnelc(x, err_id);};

  vectorize(f_fc, v, r, option);
}

void GNSTLIB::fresnelc_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_fc = 
    [](double x) {int err_id; return GNSTLIB::fresnelc(x, err_id);};

  vectorize(f_fc, n, v, r, option);
}

void GNSTLIB::fresnelc_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_fc =
    [](cmplx z) {int err_id; return GNSTLIB::fresnelc(z, err_id);};

  vectorize(f_fc, v, r, option);
}

void GNSTLIB::fresnelc_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_fc = 
    [](cmplx z) {int err_id; return GNSTLIB::fresnelc(z, err_id);};

  vectorize(f_fc, n, v, r, option);
}

void GNSTLIB::fresnels_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_fs = 
    [](double x) {int err_id; return GNSTLIB::fresnels(x, err_id);};

  vectorize(f_fs, v, r, option);
}

void GNSTLIB::fresnels_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_fs = 
    [](double x) {int err_id; return GNSTLIB::fresnels(x, err_id);};

  vectorize(f_fs, n, v, r, option);
}

void GNSTLIB::fresnels_vec(std::vector<cmplx>& v, std::vector<cmplx>& r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_fs =
    [](cmplx z) {int err_id; return GNSTLIB::fresnels(z, err_id);};

  vectorize(f_fs, v, r, option);
}

void GNSTLIB::fresnels_vec(const int n, const cmplx *v, cmplx *r, 
  unsigned short option)
{
  std::function<cmplx(cmplx)> f_fs = 
    [](cmplx z) {int err_id; return GNSTLIB::fresnels(z, err_id);};

  vectorize(f_fs, n, v, r, option);
}

void GNSTLIB::inverf_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_inverf = 
    [](double x) {int err_id; return GNSTLIB::inverf(x, err_id);};

  vectorize(f_inverf, v, r, option);
}

void GNSTLIB::inverf_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_inverf = 
    [](double x) {int err_id; return GNSTLIB::inverf(x, err_id);};

  vectorize(f_inverf, n, v, r, option);
}

void GNSTLIB::inverfc_vec(std::vector<double>& v, std::vector<double>& r, 
  unsigned short option) 
{
  std::function<double(double)> f_inverfc = 
    [](double x) {int err_id; return GNSTLIB::inverfc(x, err_id);};

  vectorize(f_inverfc, v, r, option);
}

void GNSTLIB::inverfc_vec(const int n, const double *v, double *r, 
  unsigned short option)
{
  std::function<double(double)> f_inverfc = 
    [](double x) {int err_id; return GNSTLIB::inverfc(x, err_id);};

  vectorize(f_inverfc, n, v, r, option);
}