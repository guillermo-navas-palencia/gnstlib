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
 * gnstlib_dd.hpp
 *  + includes double-double and double-complex numbers
 *  
 *  Types double-double arithmetic. This datatype consists of a pair of 64-bit
 *  IEEE double (hi, lo), where "hi" is the 64-bit floating-point value closest
 *  to the desired value (high) and "lo" is the difference between the true
 *  value and hi (lower).This type provides theoretically up to 106-bit of 
 *  precision (approximately 31 digits of accuracy)
 */

#ifndef GNSTLIB_DD_H
#define GNSTLIB_DD_H

#include <cmath>
#include <complex>

typedef std::complex<double> cmplx;

struct dd
{
  double hi, lo; // high and low order part of the number

  // constructors
  dd() { hi = 0.0; lo = 0.0; }
  dd(double h, double l) { hi = h; lo = l; }
  dd(double h) { hi = h; lo = 0.0; }
  dd(int h) { hi = (double)h; lo = 0.0; }
};

struct ddc
{
  double xhi, xlo, yhi, ylo;

  // constructors
  ddc() { xhi = 0.0; xlo = 0.0; yhi = 0.0; ylo = 0.0; }
  ddc(double xh) {xhi = xh; xlo = 0.0; yhi = 0.0; ylo = 0.0;}
  ddc(double xh, double xl, double yh, double yl) { 
      xhi = xh; xlo = xl; yhi = yh; ylo = yl; 
  }   
};

/* Utils and auxiliary functions */

// negative of double-double number (unary minus)
inline void dd_negative(dd& a) {
  a.hi = -a.hi;
  a.lo = -a.lo;
}

// negative copy of double-double number
inline void dd_negative_copy(const dd& a, dd& x) {
  x.hi = -a.hi;
  x.lo = -a.lo;
}

/* Real operations -----------------------------------------------------------*/

/* Addition inline functions */

// error-free transformatin of the sum of two floating-point numbers
// TwoSum - Knuth's algorithm (1969). This algorithm requires 6 flops.
inline void dd_add_d_d(const double a, const double b, dd& x) {
  double c;

  x.hi = a + b;
  c = x.hi - a;
  x.lo = (a - (x.hi - c)) + (b - c);
} 

// error-free transformation of the sum of two floating-point numbers.
// FastTwoSum - Dekker (1971). This algorithm requires 3 flops.
inline void dd_add_d_d_fast(const double a, const double b, dd& x) {
  x.hi = a + b;
  x.lo = (a - x.hi) + b;
}

// addition of a double-double number and a double number. This algorithm 
// requires 10 flops.
inline void dd_add_dd_d(const dd& a, const double b, dd& x) {
  dd t;

  dd_add_d_d(a.hi, b, t);
  t.lo += a.lo;
  dd_add_d_d_fast(t.hi, t.lo, x);
}

// in-place addition of a double-double number and double-number. x <- x + a
inline void dd_add_dd_d_ip(dd& x, const double a) {
  dd t;

  dd_add_d_d(x.hi, a, t);
  t.lo += x.lo;
  dd_add_d_d_fast(t.hi, t.lo, x);
}

// addition of two double-double numbers. This algorithm requires 20 flops.
inline void dd_add_dd_dd(const dd& a, const dd&b, dd& x) {
  dd s, t;

  // add two high-order parts
  dd_add_d_d(a.hi, b.hi, s);

  // add two low-order parts
  dd_add_d_d(a.lo, b.lo, t);

  // renormalize
  s.lo += t.hi;
  t.hi = s.hi + s.lo;
  s.lo -= (t.hi - s.hi);
  t.lo += s.lo;

  dd_add_d_d_fast(t.hi, t.lo, x);
}

// in-place addition of two double-double numbers. x <- x + a
inline void dd_add_dd_dd_ip(dd& x, const dd& a) {
  dd s, t;

  // add two high-order parts
  dd_add_d_d(x.hi, a.hi, s);

  // add two low-order parts
  dd_add_d_d(x.lo, a.lo, t);
  
  // renormalize
  s.lo += t.hi;
  t.hi =s.hi + s.lo;
  s.lo -= (t.hi - s.hi);
  t.lo += s.lo;

  dd_add_d_d_fast(t.hi, t.lo, x);
}

/* Multiplication inline functions */

// error-free splitting of a floating-point number into two parts
// Dekker (1971). This algorithm requires 4 flops
inline void dd_split(const double a, dd& x) {
  int sca;
  double c, frac;

  double split_max = 6.69692879491417e+299; // 2^996
  double split_val = 1.3421772900000000e+8; // 2^27 + 1

  // check exponent of a and switch to safe mode if necessary
  if (std::abs(a) > split_max) {
    // overflow after multiplication by factor
    frac = std::frexp(a, &sca);
    c = frac * split_val;
    x.hi = (frac - c) + c;
    x.lo = frac - x.hi;

    // re-scale
    x.hi = std::ldexp(x.hi, sca);
    x.lo = std::ldexp(x.lo, sca);

  } else {
    // safe split 
    c = a * split_val;
    x.hi = c - (c - a);
    x.lo = a - x.hi;
  }
}

// error-free transformation of the product of two floating-point numbers.
// TwoProduct - G. W. Veltkamp (1971). This algorithm requires 17 flops.
inline void dd_mul_d_d(const double a, const double b, dd& x) {
  dd aa, bb;

  dd_split(a, aa);
  dd_split(b, bb);

  x.hi = a * b;
  x.lo = aa.lo*bb.lo - (((x.hi - aa.hi*bb.hi) - aa.lo*bb.hi) - aa.hi*bb.lo);
}

// multiplication of double-double number by a double number. This algoritm
// requires 22 flops.
inline void dd_mul_dd_d(const dd& a, const double b, dd& x) {
  dd t;

  dd_mul_d_d(a.hi, b, t);
  t.lo += (a.lo * b);
  dd_add_d_d_fast(t.hi, t.lo, x); 
}

// in-place multiplication of a double-double number and a double number.
// x <- x * a
inline void dd_mul_dd_d_ip(dd& x, const double a) {
  dd t;

  dd_mul_d_d(x.hi, a, t);
  t.lo += (x.lo * a);
  dd_add_d_d_fast(t.hi, t.lo, x);
}

// multiplication of two double-double numbers. This algoritm requires 24 flops.
inline void dd_mul_dd_dd(const dd& a, const dd& b, dd& x) {
  dd t;

  dd_mul_d_d(a.hi, b.hi, t);
  t.lo += ((a.hi * b.lo) + (a.lo * b.hi));
  dd_add_d_d_fast(t.hi, t.lo, x);
}

// in-place multiplication of two double-double numbers. x <- x * a
inline void dd_mul_dd_dd_ip(dd& x, const dd& a) {
  dd t;

  dd_mul_d_d(x.hi, a.hi, t);
  t.lo += ((x.hi * a.lo) + (x.lo * a.hi));
  dd_add_d_d_fast(t.hi, t.lo, x);
}

// multiply double-double number by double, which must be power of 2. 
inline void dd_mul_pw2_dd(const dd& a, const double b, dd& x) {
  x.hi = a.hi * b;
  x.lo = a.lo * b;
}

/* Division and reciprocal inline functions */

// error-free transformation of the division of two double numbers. 
// This algorithm requires 30 flops.
inline void dd_div_d_d(const double a, const double b, dd& x) {
  double q1, q2;
  dd p, s;

  // first quotient approximation
  q1 = a / b;

  // compute x.hi * b
  dd_mul_d_d(q1, b, p);

  // compute a - x.hi * b
  dd_add_d_d(a, -p.hi, s);

  s.lo -= p.lo;

  // next approximation
  q2 = (s.hi + s.lo) / b;

  dd_add_d_d_fast(q1, q2,  x);
}

// division of a double-double number by double number. This algorithm requires
// 31 flops.
inline void dd_div_dd_d(const dd& a, const double b, dd& x) {
  double q1, q2;
  dd p, t;

  // first quotient approximation
  q1 = a.hi / b;

  // compute q1 * b
  dd_mul_d_d(q1, b, t);

  t.lo += a.lo;
  t.lo -= p.lo;

  // next approximation
  q2 = (t.hi + t.lo) / b;

  dd_add_d_d_fast(q1, q2, x);
}

// in-place division of double-double by double. x <- x / a.
inline void dd_div_dd_d_ip(dd& x, const double a) {
  double q1, q2;
  dd p, t;

  q1 = x.hi / a;
  dd_mul_d_d(q1, a, p);
  dd_add_d_d(x.hi, -p.hi, t);
  t.lo += x.lo;
  t.lo -= p.lo;

  q2 = (t.hi + t.lo) / a;
  dd_add_d_d_fast(q1, q2, x);
}

// reciprocal or inverse of double-double number. Use truncated Newton iteration
// based algorithm. This algorithm requires 55 flops.
inline void dd_inv_dd(const dd& a, dd& inva) {
  double x0;
  dd w, v;

  // high-order approximation
  x0 = 1.0 / a.hi;

  dd_mul_dd_d(a, x0, v);
  dd_negative(v);
  dd_add_dd_d(v, 2.0, w);
  dd_mul_dd_d(w, x0, inva);
}

// division of two double-double numbers. Compute reciprocal algorithm based on
// Newton-Raphson for 1/b and multiply by a. These two combined algorithms 
// require 55 + 24 = 79 flops. Therefore, it reduces the number of flops from 
// 100 (Bailey's QD library) to 79.
inline void dd_div_dd_dd(const dd&a, const dd& b, dd& x) {
  dd rb;

  // compute reciprocal of b
  dd_inv_dd(b, rb);

  // multiply by a
  dd_mul_dd_dd(a, rb, x);
}

// in-place division of two double-double numbers. x <- x / a
inline void dd_div_dd_dd_ip(dd& x, const dd& a) {
  dd ra;

  dd_inv_dd(a, ra);
  dd_mul_dd_dd_ip(x, ra);
}

// inverse of double number. Call dd_div_d_d
inline void dd_inv_d(const double x, dd& invx) {
  dd_div_d_d(1.0, x, invx);
}

// division of a double number by double-double. Call dd_div_dd_dd.
inline void dd_div_d_dd(const double a, const dd& b, dd& x) {
  dd aa;

  aa = a; // convert a to dd_real aa
  dd_div_dd_dd(aa, b, x);
}

/* Square */

// error-free transformation of square of double number. Faster than 
// dd_mul_d_d(a, a, x). This algorithm requires 12 flops.
inline void dd_sqr_d(const double a, dd& x) {
  dd aa;

  dd_split(a, aa);

  // square approximation
  x.hi = a * a;
    x.lo = ((aa.hi * aa.hi - x.hi) + 2.0 * aa.hi * aa.lo) + aa.lo * aa.lo;
}

// compute square of double-double number. Faster than dd_mul_dd_dd(a, a, x).
// This algorithm requires 18 flops.
inline void dd_sqr_dd(const dd& a, dd& x) {
  dd p;

  // high-order square approximation 
  dd_sqr_d(a.hi, p);
  p.lo += 2.0 * a.hi * a.lo;

  dd_add_d_d_fast(p.hi, p.lo, x);
}

/* Square root */

// square root double number. This algorithm requires 33 flops.
// Alan H. Karp and Peter Marstein (1997)
inline void dd_sqrt_d(const double a, dd& x) {
  double ax, x0, w;
  dd ax2, v;

  // high-order approximation as starting point for N-R
  x0 = 1. / std::sqrt(a);
  ax = a * x0;

  dd_sqr_d(ax, ax2);
  dd_negative(ax2);
  dd_add_dd_d(ax2, a, v);
  w = v.hi * x0 * 0.5;
  dd_add_d_d(w, ax, x);
}

// square root of double-double number. This algorithm requires 43 flops
inline void dd_sqrt_dd(const dd& a, dd& x) {
  double ax, x0, w;
  dd ax2, v;

  x0 = 1. / std::sqrt(a.hi);
  ax = a.hi * x0;
  dd_sqr_d(ax, ax2);
  dd_negative(ax2);
  dd_add_dd_dd(a, ax2, v);
  w = v.hi * x0 * 0.5;
  dd_add_d_d(w, ax, x);
}

// reciprocal of square root of double number. This algorithm requires 70 flops.
// M. Joldes et al. (2016)
inline void dd_rsqrt_d(const double a, dd& x) {
  double x0;
  dd t, v, x02;

  x0 = 1. / std::sqrt(a);
  dd_sqr_d(x0, x02);
  dd_mul_dd_d_ip(x02, a);
  dd_negative(x02);
  dd_add_dd_d(x02, 3.0, t);
  dd_mul_dd_d(t, x0, v);
  dd_mul_pw2_dd(v, 0.5, x);
}

// reciprocal of square of double-double number. This algorithm requires 77 fl.
inline void dd_rsqrt_dd(const dd& a, dd& x) {
  double x0;
  dd t, v, x02;

  x0 = 1. / std::sqrt(a.hi);
  dd_sqr_d(x0, x02);
  dd_mul_dd_dd_ip(x02, a);
  dd_negative(x02);
  dd_add_dd_d(x02, 3.0, t);
  dd_mul_dd_d(t, x0, v);
  dd_mul_pw2_dd(v, 0.5, x);
}

/* Complex operations --------------------------------------------------------*/

/* Addition inline operations */

// error-free transformation of the sum of two double complex numbers.
// TwoSumCplx - Graillat et. al (2012). This algorithm requires 12 flops.
inline void ddc_add_c_c(const cmplx a, const cmplx b, ddc& z) {
  dd s1, s2;

  // sum real and imaginary parts
  dd_add_d_d(a.real(), b.real(), s1);
  dd_add_d_d(a.imag(), b.imag(), s2);

  z.xhi = s1.hi; z.xlo = s1.lo;
  z.yhi = s2.hi; z.ylo = s2.lo;
}

// addition of double complex and double number.
inline void ddc_add_c_d(const cmplx a, const double b, ddc& z) {
  dd s1;

  // sum real part 
  dd_add_d_d(a.real(), b, s1);

  z.xhi = s1.hi; z.xlo = s1.lo;
  z.yhi = a.imag(); z.ylo = 0.0;
}

// addition of double-double complex and double complex. This algorithm requires
// 20 flops.
inline void ddc_add_ddc_c(const ddc& a, const cmplx b, ddc& z) {
  dd aa_re, aa_im, s1, s2;

  // sum real and imaginary part
  aa_re = dd(a.xhi, a.xlo);
  dd_add_dd_d(aa_re, b.real(), s1);
  aa_im = dd(a.yhi, a.ylo);
  dd_add_dd_d(aa_im, b.imag(), s2);

  z.xhi = s1.hi; z.xlo = s1.lo;
  z.yhi = s2.hi; z.ylo = s2.lo;
}

// in-place addition of double-double complex and double complex. z <- z + a
inline void ddc_add_ddc_c_ip(ddc& z, const cmplx a) {
  dd z_re, z_im, s1, s2;

  z_re = dd(z.xhi, z.xlo);
  dd_add_dd_d(z_re, a.real(), s1);
  z_im = dd(z.yhi, z.ylo);
  dd_add_dd_d(z_im, a.imag(), s2);

  z.xhi = s1.hi; z.xlo = s1.lo;
  z.yhi = s2.hi; z.ylo = s2.lo;
}

// addition of two double-double complex numbers. This algorithm requires 40 fl.
inline void ddc_add_ddc_ddc(const ddc& a, const ddc& b, ddc& z) {
  dd aa_re, aa_im, bb_re, bb_im, s1, s2;

  // sum real parts
  aa_re = dd(a.xhi, a.xlo);
  bb_re = dd(b.xhi, b.xlo);
  dd_add_dd_dd(aa_re, bb_re, s1);
  // sum imaginary parts
  aa_im = dd(a.yhi, a.ylo);
  bb_im = dd(b.yhi, b.ylo);
  dd_add_dd_dd(aa_im, bb_im, s2);

  z.xhi = s1.hi; z.xlo = s1.lo;
  z.yhi = s2.hi; z.ylo = s2.lo; 
}

// in-place addition of two double-double complex numbers. z <- z + a
inline void ddc_add_ddc_ddc_ip(ddc& z, const ddc& a) {
  dd z_re, z_im, aa_re, aa_im, s1, s2;

  z_re = dd(z.xhi, z.xlo);
  aa_re = dd(a.xhi, a.xlo);
  dd_add_dd_dd(z_re, aa_re, s1);
  z_im = dd(z.yhi, z.ylo);
  aa_im = dd(a.yhi, a.ylo);
  dd_add_dd_dd(z_im, aa_im, s2);

  z.xhi = s1.hi; z.xlo = s1.lo;
  z.yhi = s2.hi; z.ylo = s2.lo;   
}


/* Complex multiplication */

// error-free transformation of the product of two complex numbers. 
// TwoProductCplxSingeSplitting - Graillat et al. (2012). This algorithm 
// requires 92 flops.

inline void ddc_mul_c_c(const cmplx x, const cmplx y, ddc& z) {
  double a, b, c, d;
  dd a_, b_, c_, d_, ac, bd , bc, ad, xx, yy;

  a =  x.real(); b = x.imag(); c = y.real(); d = y.imag();

  // split a, b, c, d
  dd_split(a, a_);
  dd_split(b, b_);
  dd_split(c, c_);
  dd_split(d, d_);

  // compute a*c
  ac.hi = a * c;
  ac.lo = a_.lo*c_.lo - (((ac.hi-a_.hi*c_.hi)-a_.lo*c_.hi) - a_.hi*c_.lo);

  // compute b * d
  bd.hi = b * d;
  bd.lo = b_.lo*d_.lo - (((bd.hi-b_.hi*d_.hi)-b_.lo*d_.hi) - b_.hi*d_.lo);

  // compute b * c
  bc.hi = b * c;
  bc.lo = b_.lo*c_.lo - (((bc.hi-b_.hi*c_.hi)-b_.lo*c_.hi) - b_.hi*c_.lo);

  // compute a * d
  ad.hi = a * d;
  ad.lo = a_.lo*d_.lo - (((ad.hi-a_.hi*d_.hi)-a_.lo*d_.hi) - a_.hi*d_.lo);

  dd_add_d_d(ac.hi, -bd.hi, xx);
  dd_add_d_d(bc.hi, ad.hi, yy);

  z.xhi = xx.hi; z.xlo = xx.lo + ac.lo - bd.lo;
  z.yhi = yy.hi; z.ylo = yy.lo + ad.lo + bc.lo;
}

// multiplication of double-double complex by a complex number. This algorithm
// requires 112 flops.
inline void ddc_mul_ddc_c(const ddc& x, const cmplx y, ddc& z) {
  double c, d;
  dd aa, aahi, bb, bbhi, u1, u2, u3, u4, ac, bd, bc, ad, cc, dd_, xx, yy;

  aa = dd(x.xhi, x.xlo);
  bb = dd(x.yhi, x.ylo);
  c = y.real();
  d = y.imag();

  // split c and d
  dd_split(aa.hi, aahi);
  dd_split(bb.hi, bbhi);
  dd_split(c, cc);
  dd_split(d, dd_);

  // compute a * c
  u1.hi = aa.hi * c;
  u1.lo = u1.lo*cc.lo - (((u1.hi-aahi.hi*cc.hi)-aahi.lo*cc.hi)-aahi.hi*cc.lo);
  u1.lo += aa.lo * c; 
  dd_add_d_d_fast(u1.hi, u1.lo, ac);

  // compute b * d
  u2.hi = bb.hi * d;
  u2.lo = bbhi.lo*dd_.lo - (((u2.hi-bbhi.hi*dd_.hi)-bbhi.lo*dd_.hi)
    -bbhi.hi*dd_.lo);
  u2.lo += bb.lo * d;
  dd_add_d_d_fast(u2.hi, u2.lo, bd);

  // compute b * c
  u3.hi = bb.hi * c;
  u3.lo = bbhi.lo*cc.lo - (((u3.hi-bbhi.hi*cc.hi)-bbhi.lo*cc.hi)
    -bbhi.hi*cc.lo);
  u3.lo += bb.lo * c;
  dd_add_d_d_fast(u3.hi, u3.lo, bc);

  // compute a * d
  u4.hi = aa.hi * d;
  u4.lo = aahi.lo*dd_.lo - (((u4.hi-aahi.hi*dd_.hi)-aahi.lo*dd_.hi)
    -aahi.hi*dd_.lo);
  u4.lo += aa.lo * d;
  dd_add_d_d_fast(u4.hi,u4.lo, ad);

  dd_negative(bd);
  dd_add_dd_dd(ac, bd, xx);
  dd_add_dd_dd(bc, ad, yy);

  z.xhi = xx.hi; z.xlo = xx.lo;
  z.yhi = yy.hi; z.ylo = yy.lo;
}

// inline void ddc_mul_ddc_c_ip(const cmplx a, ddc& z) {}

// inline void ddc_mul_ddc_ddc(const ddc& x, const ddc& y, ddc& z) {}

// inline void ddc_mul_ddc_ddc_ip(const ddc& a, ddc& z) {}

/* Complex Division */
#endif /* GNSTLIB_DD_H */