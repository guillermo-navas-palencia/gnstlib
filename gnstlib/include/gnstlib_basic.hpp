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
 *  + includes elementary functions
 *    - exponential and logarithmic functions
 *    - trigonometric functions
 *    - hyperbolic functions
 *    - Power functions
 */

#ifndef GNSTLIB_BASIC_H
#define GNSTLIB_BASIC_H

#include <complex>

namespace GNSTLIB 
{

/* Exponential and logarithmic functions */
/******************************************************************************/
// compute exp(x)
double exp(const double x);
// compute exp2(x) = 2^x
double exp2(const double x);
// compute exp10(x) = 10^x
double exp10(const double x);
// compute exp(x) - 1
double expm1(const double x);
// compute (exp(x) - 1) / x
double expm1x(const double x);
// compute (exp(x) - 1 - x) / (0.5 * x * x)
double expm1mx(const double x);

// compute log(x)
double log(const double x);
// compute log2(x)
double log2(const double x);
// compute log10(x)
double log10(const double x);
// compute log(1 + x) = log1p(x)
double log1p(const double x);
// compute log(1 + x) - x = log1p(x) - x
double log1pmx(const double x); 
// compute log(1 - exp(x))
double log1mexp(const double x, int& err_id);
// compute log(1 + exp(x))
double log1pexp(const double x);

//compute exp(z)
std::complex<double> exp(const std::complex<double> z);
//compute log(z)
std::complex<double> log(const std::complex<double> z);

/* Trigonometric functions */
/******************************************************************************/
// compute cos(x)
double cos(const double x);
// compute sin(x)
double sin(const double x);
// compute tan(x)
double tan(const double x);
// compute sec(x)
double sec(const double x);
// compute csc(x)
double csc(const double x);
// compute cot(x)
double cot(const double x);

// compute acos(x)
double acos(const double x);
// compute asin(x)
double asin(const double x);
// compute atan(x)
double atan(const double x);
// còmpute atan2(y, x)
double atan2(const double y, const double x);
// compute asec(x)
double asec(const double x);
// compute acsc(x)
double acsc(const double x);
// compute acot(x)
double acot(const double x);

// compute cos(pi * x)
double cospi(const double x);
// compute sin(pi * x)
double sinpi(const double x);

// compute cos(z)
std::complex<double> cos(const std::complex<double> z);
// compute sin(z)
std::complex<double> sin(const std::complex<double> z);
// compute tan(z)
std::complex<double> tan(const std::complex<double> z);
// compute sec(z)
std::complex<double> sec(const std::complex<double> z);
// compute csc(z)
std::complex<double> csc(const std::complex<double> z);
// compute cot(z)
std::complex<double> cot(const std::complex<double> z);

// compute acos(z)
std::complex<double> acos(const std::complex<double> z);
// compute asin(z)
std::complex<double> asin(const std::complex<double> z);
// compute atan(z)
std::complex<double> atan(const std::complex<double> z);
// compute asec(z)
std::complex<double> asec(const std::complex<double> z);
// compute acsc(z)
std::complex<double> acsc(const std::complex<double> z);
// compute acot(z)
std::complex<double> acot(const std::complex<double> z);

// compute cos(pi * z)
std::complex<double> cospi(const std::complex<double> z);
// compute sin(pi * z)
std::complex<double> sinpi(const std::complex<double> z);

/* Hyperbolic functions */
/******************************************************************************/
// compute cosh(x)
double cosh(const double x);
// compute sinh(x)
double sinh(const double x);
// compute tanh(x)
double tanh(const double x);
// compute sech(x)
double sech(const double x);
// compute csch(x)
double csch(const double x);
// compute coth(x)
double coth(const double x);

// compute acosh(x)
double acosh(const double x);
// compute asinh(x)
double asinh(const double x);
// compute atanh(x)
double atanh(const double x);
// compute asech(x)
double asech(const double x);
// compute acsch(x)
double acsch(const double x);
// compute acoth(x)
double acoth(const double x);

// compute cosh(z)
std::complex<double> cosh(const std::complex<double> z);
// compute sinh(z)
std::complex<double> sinh(const std::complex<double> z);
// compute tanh(z)
std::complex<double> tanh(const std::complex<double> z);
// compute sech(z)
std::complex<double> sech(const std::complex<double> z);
// compute csch(z)
std::complex<double> csch(const std::complex<double> z);
// compute coth(z)
std::complex<double> coth(const std::complex<double> z);

// compute acosh(z)
std::complex<double> acosh(const std::complex<double> z);
// compute asinh(z)
std::complex<double> asinh(const std::complex<double> z);
// compute atanh(z)
std::complex<double> atanh(const std::complex<double> z);
// compute asech(z)
std::complex<double> asech(const std::complex<double> z);
// compute acsch(z)
std::complex<double> acsch(const std::complex<double> z);
// compute acoth(z)
std::complex<double> acoth(const std::complex<double> z);

/* Power functions */
/******************************************************************************/
// compute pow(x, y)
double pow(const double x, const double y);
// compute sqrt(x)
double sqrt(const double x);
// compute hypot(x, y)
double hypot(const double x, const double y);
// compute cbrt(x)
double cbrt(const double x);
// compute powm1(x, y) = x^y - 1
double powm1(const double x, const double y, int& err_id);

// compute power(cx, cy)
std::complex<double> pow(const std::complex<double> x, 
                         const std::complex<double> y);
// compute sqrt(z)
std::complex<double> sqrt(const std::complex<double> z);

} // namespace GNSTLIB

#endif /* GNSTLIB_BASIC_H */