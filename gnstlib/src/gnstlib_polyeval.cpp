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
#include <vector>
#include <complex>

#include "gnstlib_dd.hpp"

typedef std::complex<double> cmplx;

/**
 * Polynomial evaluation, real coefficients and complex variable
 */
cmplx polyeval(const std::vector<double>& v, const cmplx z) {
  unsigned int j, n;
  double a = v.back();
  double b = v.end()[-1];
  double r = 2.0 * real(z);
  double s = real(z) * real(z) + imag(z) * imag(z);
  double tmp;

  n = static_cast<int>(v.size());
  for (j = 2; j <= n; j++) {
    tmp = b;
    b = v[n-j] - s * a;
    a = r * a + tmp;
  } 
  return z * a + b;
}

/**
 * Polynomial evaluation, real coefficients and complex variable 
 * (accurate version)
 */
ddc polyeval_acc(const std::vector<double>& v, const cmplx z) {
  unsigned int j, n;
  dd a = v.back();
  dd b = v.end()[-1];
  double r = 2.0 * real(z);
  double s = real(z) * real(z) + imag(z) * imag(z);
  dd tmp;
  ddc az, azz, bz, rz;

  n = static_cast<int>(v.size());
  for (j = 2; j <= n; j++) {
    tmp = b;
    dd_mul_dd_d(a, s, b);
    dd_negative(b);
    dd_add_dd_d_ip(b, v[n-j]);
    dd_mul_dd_d_ip(a, r);
    dd_add_dd_dd_ip(a, tmp);
  }
  az = ddc(a.hi, a.lo, 0.0, 0.0);
  bz = ddc(b.hi, b.lo, 0.0, 0.0);
  ddc_mul_ddc_c(az, z, azz);
  ddc_add_ddc_ddc(azz, bz, rz);
  return rz;
}