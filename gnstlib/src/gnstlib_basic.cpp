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
#include "gnstlib_dd.hpp"
#include "gnstlib_utils.hpp"

typedef std::complex<double> cmplx;

namespace GNSTLIB
{
/* Exponential and logarithmic functions */

// compute exp(x)
double exp(const double x) 
{
  double ans = std::exp(x);
  clean_result(ans);
  return ans;
}

// compute exp(z)
cmplx exp(const cmplx z) 
{
  cmplx ans = std::exp(z);
  clean_result_complex(ans);
  return ans;
}

// compute exp2(x): 2^x
double exp2(const double x)
{
  double ans = std::exp2(x);
  clean_result(ans);
  return ans;
}

// compute exp10(x): 10^x
// based on musl - an implementation of the standard library for Linux-based
double exp10(const double x)
{
  const double maxexp10 = 308.2547155599167439;

  // special cases
  if (std::isnan(x))
    return nan;
  if (x > maxexp10)
    return inf;
  if (x < -maxexp10)
    return 0.0;

  const std::vector<double> p10 = {1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 
    1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 
    1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15};

  double n, y = std::modf(x, &n);
  union {double f; uint64_t i;} u = {n};
  if ((u.i >> 52 & 0x7ff) < 0x3ff+4)
  {
    if (!y)
      return p10[static_cast<int>(n) + 15];
    y = std::exp2(3.32192809488736234787031942948939 * y);
    return y * p10[static_cast<int>(n) + 15];
  }
  return std::pow(10.0, x);
}

// compute exp(x) - 1
double expm1(const double x)
{
  double ans = std::expm1(x);
  clean_result(ans);
  return ans;
}

// compute (exp(x) - 1) / x
double expm1x(const double x) 
{
  double t;

  if (x == 0.0) 
    return 1.0;
  else if ((x < -0.69) || (x > 0.4))
    return std::expm1(x) / x;
  else {
    t = 0.5 * x;
    return std::exp(t) * std::sinh(t) / t;
  }
}

// compute (exp(x) - 1 - x) / (0.5 * x * x)
double expm1mx(const double x) 
{
  double t, t2;

  if (x == 0.0)
    return 1.0;
  else if (std::abs(x) > 0.9) 
    return 2.0 * (std::expm1(x) - x) / (x * x);
  else {
    t = std::sinh(0.5 * x);
    t2 = t * t;
    return 2.0 * (2.0 * t2 + (2.0 * t * std::sqrt(1.0 + t2) - x)) / (x * x);
  }
}

// compute log(x)
double log(const double x) 
{
  double ans = std::log(x);
  clean_result(ans);
  return ans;
}

// compute log(z)
cmplx log(const cmplx z) 
{
  cmplx ans = std::log(z);
  clean_result_complex(ans);
  return ans;
}

// compute log2(x)
double log2(const double x) 
{
  double ans = std::log2(x);
  clean_result(ans);
  return ans;
}

// compute log10(x)
double log10(const double x) 
{
  double ans = std::log10(x);
  clean_result(ans);
  return ans;
}

// compute log1p(x): log(1+x)
double log1p(const double x)
{
  double ans = std::log1p(x);
  clean_result(ans);
  return ans; 
}

// compute log(1 + x) - x = log1p(x) - x
double log1pmx(const double x) 
{
  double e2, r, s, y0, z;

  z = std::log1p(x);
  y0= z - x;
  e2 = expm1mx(z);
  s = 0.5 * e2 * z * z;
  r = (s + y0) / (s + 1.0 + z);
  return y0 - r * (6.0 - r) / (6.0 - 4.0 * r);
}

// compute log(1 + exp(x))
// https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
double log1pexp(const double x) 
{
  if (x <= -37.0)
    return std::exp(x);
  else if (x > -37 && x <= 18.0)
    return std::log1p(std::exp(x));
  else if (x > 18.0 && x <= 33.3)
    return x + std::exp(-x);
  else
    return x;
}

// compute log(1 - exp(x)), x < 0.0
// https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
double log1mexp(const double x, int& err_id) 
{
  err_id = 0;

  if (x == 0.0)
    return -inf;

  if (x >= 0.0) {
    err_id = 1;
    return nan;
  }
  
  double a = -x;
  if (a <= constants::log_2 && a > 0.0)
    return std::log(-std::expm1(x));
  else
    return std::log1p(-std::exp(x));
}

/* Trigonometric functions */

// compute sin(x)
double sin(const double x)
{
  double ans = std::sin(x);
  clean_result(ans);
  return ans;
}

// compute sin(z)
cmplx sin(const cmplx z)
{
  cmplx ans = std::sin(z);
  clean_result_complex(ans);
  return ans;
}

// compute cos(x)
double cos(const double x)
{
  double ans = std::cos(x);
  clean_result(ans);
  return ans;
}

// compute cos(z)
cmplx cos(const cmplx z)
{
  cmplx ans = std::cos(z);
  clean_result_complex(ans);
  return ans;
}

// compute tan(x)
double tan(const double x)
{
  double ans = std::tan(x);
  clean_result(ans);
  return ans;
}

// compute tan(z)
cmplx tan(const cmplx z)
{
  cmplx ans = std::tan(z);
  clean_result_complex(ans);
  return ans;
}

// compute sec(x): sec(x) = 1 / cos(x)
double sec(const double x)
{
  double ans = 1.0 / std::cos(x);
  clean_result(ans);
  return ans;
}

// compute sec(z): sec(z) = 1 / cos(z)
cmplx sec(const cmplx z)
{
  cmplx ans = 1.0 / std::cos(z);
  clean_result_complex(ans);
  return ans;
}

// compute csc(x): csc(x) = 1 / sin(x)
double csc(const double x)
{
  double ans = 1.0 / std::sin(x);
  clean_result(ans);
  return ans;
}

// compute csc(z): csc(z) = 1 / sin(z)
cmplx csc(const cmplx z)
{
  cmplx ans = 1.0 / std::sin(z);
  clean_result_complex(ans);
  return ans;
}

// compute cot(x): 1 / tan(x) = cos(x) / sin(x)
double cot(const double x)
{
  double ans = 1.0 / std::tan(x);
  clean_result(ans);
  return ans; 
}

// compute cot(z): 1 / tan(z) = cos(z) / sin(z)
cmplx cot(const cmplx z)
{
  cmplx ans = 1.0 / std::tan(z);
  clean_result_complex(ans);
  return ans; 
}

// compute acos(x)
double acos(const double x)
{
  double ans = std::acos(x);
  clean_result(ans);
  return ans;
}

// compute acos(z)
cmplx acos(const cmplx z)
{
  cmplx ans = std::acos(z);
  clean_result_complex(ans);
  return ans;
}

// compute asin(x)
double asin(const double x)
{
  double ans = std::asin(x);
  clean_result(ans);
  return ans;
}

// compute asin(z)
cmplx asin(const cmplx z)
{
  cmplx ans = std::asin(z);
  clean_result_complex(ans);
  return ans;
}

// compute atan(x)
double atan(const double x)
{
  double ans = std::atan(x);
  clean_result(ans);
  return ans;
}

// compute atan(z)
cmplx atan(const cmplx z)
{
  cmplx ans = std::atan(z);
  clean_result_complex(ans);
  return ans;
}

// compute atan2(x, y): signed angle between the positive x-axis and point (x,y)
double atan2(const double y, const double x)
{
  double ans = std::atan2(y, x);
  clean_result(ans);
  return ans;
}

// compute asec(x)
double asec(const double x)
{
  double ans = std::acos(1.0 / x);
  clean_result(ans);
  return ans;
}

// compute asec(z)
cmplx asec(const cmplx z)
{
  cmplx ans = std::acos(1.0 / z);
  clean_result_complex(ans);
  return ans;
} 

// compute acsc(x)
double acsc(const double x)
{
  double ans = std::asin(1.0 / x);
  clean_result(ans);
  return ans;
}

// compute acsc(z)
cmplx acsc(const cmplx z)
{
  cmplx ans = std::asin(1.0 / z);
  clean_result_complex(ans);
  return ans;
}

// compute acot(x)
double acot(const double x)
{
  double ans = std::atan(1.0 / x);
  clean_result(ans);
  return ans;
}

// compute acot(z)
cmplx acot(const cmplx z)
{
  cmplx ans = std::atan(1.0 / z);
  clean_result_complex(ans);
  return ans;
}

// compute sin(pi * x) 
// based on boost implementation
double sinpi(const double x)
{
  double ans, rem;
  double u = x;
  bool invert = false;

  if (x < 0.0)
    return -sinpi(-x);

  if (x == 0.5)
    return 1.0;

  if (x < 0.5)
    return std::sin(constants::pi * x);

  if (x < 1.0) {
    invert = true;
    u = -u;
  } else 
    invert = false;

  rem = std::floor(u);
  if (int(std::trunc(rem)) & 1)
    invert = !invert;

  rem = u - rem;
  if (rem > 0.5)
    rem = 1.0 - rem;

  ans = std::sin(constants::pi * rem);

  return invert ? -ans : ans;
}

// compute cos(pi * x)
// based on boost implementation
double cospi(const double x) 
{
  double ans, rem;
  double u = x;
  bool invert = false;

  if (x == 0.5)
    return 0.0;

  if (std::abs(x) < 0.25)
    return std::cos(constants::pi * x);

  if (x < 0.0)
    u = -u;

  rem = std::floor(u);
  if (int(std::trunc(rem)) & 1)
    invert = !invert;
  
  rem = u - rem;
  if (rem > 0.5) {
    rem = 1.0 - rem;
    invert = !invert;
  }
  
  if (rem > 0.25) {
    rem = 0.5 - rem;
    ans = sinpi(rem);
  } else {
    ans = std::cos(constants::pi * rem);
  }

  return invert ? -ans : ans;
}

// compute sin(pi * z)
// based on scipy implementation scipy/special/_trig.pxd
cmplx sinpi(const cmplx z) 
{
  double x = real(z);
  double piy = constants::pi * imag(z);
  double abspiy = std::abs(piy);
  double sinpix = sinpi(x);
  double cospix = cospi(x);
  double exphpiy, coshfac, sinhfac;
  cmplx ans;

  if (abspiy < 700.0) {
    ans = cmplx(sinpix * std::cosh(piy), cospix * std::sinh(piy));
    if (x == std::floor(x)) // i * sinh(imag(z) * pi)
      return cmplx(0.0, imag(ans));
    else 
      return ans;
  }

  exphpiy = std::exp(abspiy/2);
  if (std::isinf(exphpiy)) {
    if (sinpix == 0.0)
      coshfac = std::copysign(0.0, sinpix);
    else
      coshfac = std::copysign(inf, sinpix);
    if (cospix == 0.0)
      sinhfac = std::copysign(0.0, cospix);
    else
      sinhfac = std::copysign(inf, cospix);
    return cmplx(coshfac, sinhfac);
  }

  coshfac = 0.5 * sinpix * exphpiy;
  sinhfac = 0.5 * cospix * exphpiy;
  
  ans = cmplx(coshfac * exphpiy, sinhfac * exphpiy);
  if (x == std::floor(x)) // i * sinh(imag(z) * pi)
    return cmplx(0.0, imag(ans));
  else 
    return ans;
}

// compute cos(pi * z)
// based on scipy implementation scipy/special/_trig.pxd
cmplx cospi(const cmplx z) 
{
  double x = real(z);
  double piy = constants::pi * imag(z);
  double abspiy = std::abs(piy);
  double sinpix = sinpi(x);
  double cospix = cospi(x);
  double exphpiy, coshfac, sinhfac;
  cmplx ans;

  if (abspiy < 700.0) {
    ans = cmplx(cospix * std::cosh(piy), -sinpix * std::sinh(piy));
    if (x == std::floor(x)) // -cosh(imag(z) * pi)
      return cmplx(real(ans), 0.0);
    else
      return ans;
  }

  exphpiy = std::exp(abspiy/2);
  if (std::isinf(exphpiy)) {
    if (sinpix == 0.0)
      coshfac = std::copysign(0.0, sinpix);
    else
      coshfac = std::copysign(inf, sinpix);
    if (cospix == 0.0)
      sinhfac = std::copysign(0.0, cospix);
    else
      sinhfac = std::copysign(inf, cospix);
    return cmplx(coshfac, sinhfac);
  }

  coshfac = 0.5 * cospix * exphpiy;
  sinhfac = 0.5 * sinpix * exphpiy;
  ans = cmplx(coshfac * exphpiy, sinhfac * exphpiy);
  if (x == std::floor(x)) // -cosh(imag(z) * pi)
    return cmplx(real(ans), 0.0);
  else
    return ans;
}

/* Hyperbolic functions */

// compute cosh(x)
double cosh(const double x)
{
  double ans = std::cosh(x);
  clean_result(ans);
  return ans;
}

// compute cosh(z)
cmplx cosh(const cmplx z)
{
  cmplx ans = std::cosh(z);
  clean_result_complex(ans);
  return ans;
}

// compute sinh(x)
double sinh(const double x)
{
  double ans = std::sinh(x);
  clean_result(ans);
  return ans;
}

// compute sinh(z)
cmplx sinh(const cmplx z)
{
  cmplx ans = std::sinh(z);
  clean_result_complex(ans);
  return ans;
}

// compute tanh(x)
double tanh(const double x)
{
  double ans = std::tanh(x);
  clean_result(ans);
  return ans;
}

// compute tanh(z)
cmplx tanh(const cmplx z)
{
  cmplx ans = std::tanh(z);
  clean_result_complex(ans);
  return ans;
}

// compute sech(z)
double sech(const double x)
{
  double ans = 1.0 / std::cosh(x);
  clean_result(ans);
  return ans;
}

// compute sech(z)
cmplx sech(const cmplx z)
{
  cmplx ans = 1.0 / std::cosh(z);
  clean_result_complex(ans);
  return ans;
}

// compute csch(x)
double csch(const double x)
{
  double ans = 1.0 / std::sinh(x);
  clean_result(ans);
  return ans;
}

// compute csch(z)
cmplx csch(const cmplx z)
{
  cmplx ans = 1.0 / std::sinh(z);
  clean_result_complex(ans);
  return ans;
}

// compute coth(x)
double coth(const double x)
{
  double ans = 1.0 / std::tanh(x);
  clean_result(ans);
  return ans;
}

// compute coth(z)
cmplx coth(const cmplx z)
{
  cmplx ans = 1.0 / std::tanh(z);
  clean_result_complex(ans);
  return ans;
}

// compute acosh(x)
double acosh(const double x)
{
  double ans = std::acosh(x);
  clean_result(ans);
  return ans;
}

// compute acosh(z)
cmplx acosh(const cmplx z)
{
  cmplx ans = std::acosh(z);
  clean_result_complex(ans);
  return ans;
}

// compute asinh(x)
double asinh(const double x)
{
  double ans = std::asinh(x);
  clean_result(ans);
  return ans;
}

// compute asinh(z)
cmplx asinh(const cmplx z)
{
  cmplx ans = std::asinh(z);
  clean_result_complex(ans);
  return ans;
}

// compute atanh(x)
double atanh(const double x)
{
  double ans = std::atanh(x);
  clean_result(ans);
  return ans;
}

// compute atanh(z)
cmplx atanh(const cmplx z)
{
  cmplx ans = std::atanh(z);
  clean_result_complex(ans);
  return ans;
}

// compute asech(x)
double asech(const double x)
{
  double ans = std::acosh(1.0 / x);
  clean_result(ans);
  return ans;
}

// compute asech(z)
cmplx asech(const cmplx z)
{
  cmplx ans = std::acosh(1.0 / z);
  clean_result_complex(ans);
  return ans;
}

// compute acsch(x)
double acsch(const double x)
{
  double ans = std::asinh(1.0 / x);
  clean_result(ans);
  return ans;
}

// compute acsch(z)
cmplx acsch(const cmplx z)
{
  cmplx ans = std::asinh(1.0 / z);
  clean_result_complex(ans);
  return ans;
}

// compute acoth(x)
double acoth(const double x)
{
  double ans = std::atanh(1.0 / x);
  clean_result(ans);
  return ans;
}

// compute acoth(z)
cmplx acoth(const cmplx z)
{
  cmplx ans = std::atanh(1.0 / z);
  clean_result_complex(ans);
  return ans;
}


/* Power functions */

// compute pow(x, y)
double pow(const double x, const double y)
{
  double ans = std::pow(x, y);
  clean_result(ans);
  return ans;   
}

cmplx pow(const cmplx x, const cmplx y)
{
  cmplx ans = std::pow(x, y);
  clean_result_complex(ans);
  return ans;
}

// compute sqrt(x)
double sqrt(const double x)
{
  double ans = std::sqrt(x);
  clean_result(ans);
  return ans;   
}

// compute sqrt(z)
cmplx sqrt(const cmplx z)
{
  cmplx ans = std::sqrt(z);
  clean_result_complex(ans);
  return ans;
}

// compute hypot(x, y)
double hypot(const double x, const double y)
{
  double ans = std::hypot(x, y);
  clean_result(ans);
  return ans;
}

// compute cbrt(x): cubic root
double cbrt(const double x)
{
  double ans = std::cbrt(x);
  clean_result(ans);
  return ans;
}

// compute x^y - 1: useful when y is very small or x is close to 1
double powm1(const double x, const double y, int& err_id)
{
  dd ans;

  err_id = 0;

  // special cases
  if (x < 0.0 && std::trunc(y) != y) {
    err_id = 1;
    return nan;
  }

  if (x > 0.0) {
    if ((std::abs(y * (x - 1.0)) < 0.5) || std::abs(y) < 0.2) {
      double l = y * std::log(x);
      if (l < 0.5)
        return std::expm1(l);
      else 
        return inf; // overflow
    }
  } else {
    if (std::trunc(y / 2) == y / 2)
      return powm1(-x, y, err_id);
  }
  // accurate addition of two double numbers
  dd_add_d_d(std::pow(x, y), -1.0, ans);
  return ans.hi;
}

} // namespace GNSTLIB