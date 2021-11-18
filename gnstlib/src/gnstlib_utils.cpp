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
#include <complex>

#include "gnstlib.hpp"

typedef std::complex<double> cmplx;

/* FP error handler */
int gnstlib_fp_error_handler(double& result)
{
  int err_id;

  switch(std::fpclassify(result)) {
    case FP_NORMAL:    err_id = 0; break;
    case FP_INFINITE:  err_id = 1; break;
    case FP_ZERO:      err_id = 2; break;
    case FP_SUBNORMAL: err_id = 3; break;
    case FP_NAN:       err_id = -1; break;
    default:           err_id = -1; break;
  }

  if (err_id == 2)
    result = std::abs(result);
  else if (err_id == -1)
    result = nan;

  return err_id;
}

int gnstlib_fp_error_handler(cmplx& result)
{
  double re_result = std::real(result);
  double im_result = std::imag(result);

  int err_id, err_id_re = 0, err_id_im = 0;

  if (re_result)
    err_id_re = gnstlib_fp_error_handler(re_result);
  if (im_result)
    err_id_im = gnstlib_fp_error_handler(im_result);

  // infinity
  if (err_id_re == 2 || err_id_im == 2)
    err_id = 2;
  else
    err_id = std::max(err_id_re, err_id_im);

  result = cmplx(re_result, im_result);
  return err_id;
}

/* Output utils */
void clean_result(double& result) 
{
  // detect subnormal/denormal numbers and return 0 to avoid loss of precision
  if (std::fpclassify(result) == FP_SUBNORMAL)
    result = 0.0;
}

void clean_result_complex(cmplx& result)
{
  // detect subnormal/denormal numbers in real and imaginary part and return 0
  // to avoid returning results with significant loss of precision
  double _re, _im;
  _re = real(result); _im = imag(result);

  if (std::fpclassify(real(result)) == FP_SUBNORMAL)
    _re = 0.0;

  if (std::fpclassify(imag(result)) == FP_SUBNORMAL)
    _im = 0.0;

  result = cmplx(_re, _im);
}

// stopping criterion for-loops involving complex numbers
bool convergence(const cmplx s, const cmplx sp) 
{
  const double tol = 0.1 * epsilon;
  double re_a, im_a, re_b, im_b;

  re_a = real(s);
  im_a = imag(s);
  re_b = real(sp);
  im_b = imag(sp);

  return std::abs(re_a - re_b) <= tol * std::abs(re_a) &&
  std::abs(im_a - im_b) <= tol * std::abs(im_a);
}

// stopping criterion for-loops involving real numbers
bool convergence(const double s, const double sp)
{
  return std::abs((s - sp) / s) < 0.5 * epsilon;
}

/* Other mathematical functions */

// return sign
int sign(const double x) 
{
  return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}