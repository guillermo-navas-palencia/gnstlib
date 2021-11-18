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
 * gnstlib_basic.hpp
 *  + includes polynomial evaluation algorithm
 */

#ifndef GNSTLIB_POLYEVAL_H
#define GNSTLIB_POLYEVAL_H 

#include <complex>
#include <vector>

#include "gnstlib_dd.hpp"

/**
 * Chebyshev summation - Clenshaw algorithm
 */
inline double chepolsum(const std::vector<double>& v, const double x)
{
  int k;
  double u0, u1, u2, s, tt;

  int n = static_cast<int>(v.size());
  u0 = u1 = u2 = 0.0;
  tt = x + x;

  for (k = n; k --> 0; ) {
    u2 = u1;
    u1 = u0;
    u0 = tt * u1 - u2 + v[k];
  }

  s = 0.5 * (u0 - u2);
  return s;
}


/**
 * Polynomial evaluation - horner rule
 */
inline double polyeval(const std::vector<double>& v, const double x)
{
  int i;
  int n = static_cast<int>(v.size());
  double r = 0.0;

  for (i = n - 1; i >= 0; i--)
    r = r * x + v[i];

  return r;
}

// C array compatibility
inline double polyeval(const double *v, const int n, const double x) 
{
  int i;
  double r = 0.0;

  for (i = n - 1; i >= 0; i--)
    r = r * x + v[i];

  return r;
}

/**
 * Polynomial evaluation using double-double numbers (accurate version)
 */ 
inline dd polyeval_acc(const double *v, const int n, const double x) 
{
  int i;
  dd r = 0.0;
  dd u;

  for (i = n - 1; i >= 0; i--) {
    dd_mul_dd_d_ip(r, x);
    dd_add_dd_d_ip(r, v[i]);
  }

  return r;
}

inline double polyeval_acc(const std::vector<double>& v, const double x, dd& r);

std::complex<double> polyeval(const std::vector<double>& v, 
                              const std::complex<double> z);
std::complex<double> polyeval_acc(const std::vector<double>& v, 
                                  const std::complex<double> z);

// LOOP-UNROLLING polynomial evaluation
inline double polyeval_unroll3(const double *v, const double x)
{
  return (v[2] * x + v[1]) * x + v[0];
}

inline double polyeval_unroll4(const double *v, const double x)
{
  return ((v[3] * x + v[2]) * x + v[1]) * x + v[0];
}

inline double polyeval_unroll5(const double *v, const double x)
{
  return (((v[4] * x + v[3]) * x + v[2]) * x + v[1]) * x + v[0];
}

inline double polyeval_unroll6(const double *v, const double x)
{
  return ((((v[5] * x + v[4]) * x + v[3]) * x + v[2]) * x + v[1]) * x + v[0];
}

inline double polyeval_unroll7(const double *v, const double x)
{
  return (((((v[6] * x + v[5]) * x + v[4]) * x + v[3]) * x + v[2]) * x + 
    v[1]) * x + v[0];
}

inline double polyeval_unroll8(const double *v, const double x)
{
  return ((((((v[7] * x + v[6]) * x + v[5]) * x + v[4]) * x + v[3]) * x + 
    v[2]) * x + v[1]) * x + v[0];
}

inline double polyeval_unroll9(const double *v, const double x)
{
  return (((((((v[8] * x + v[7]) * x + v[6]) * x + v[5]) * x + v[4]) * x + 
    v[3]) * x + v[2]) * x + v[1]) * x + v[0];
}

inline double polyeval_unroll10(const double *v, const double x)
{
  return ((((((((v[9] * x + v[8]) * x + v[7]) * x + v[6]) * x + v[5]) * x + 
    v[4]) * x + v[3]) * x + v[2]) * x + v[1]) * x + v[0];
}

inline double polyeval_unroll11(const double *v, const double x)
{
  return (((((((((v[10] * x + v[9]) * x + v[8]) * x + v[7]) * x + v[6]) * x + 
    v[5]) * x + v[4]) * x + v[3]) * x + v[2]) * x + v[1]) * x + v[0];
}

inline double polyeval_unroll12(const double *v, const double x)
{
  return ((((((((((v[11] * x + v[10]) * x + v[9]) * x + v[8]) * x + v[7]) * x + 
    v[6]) * x + v[5]) * x + v[4]) * x + v[3]) * x + v[2]) * x + v[1]) * x + v[0];
}

#endif /* GNSTLIB_POLYEVAL_H */