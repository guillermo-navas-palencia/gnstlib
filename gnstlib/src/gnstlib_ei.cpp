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
#include <omp.h>

#include "gnstlib.hpp"
#include "gnstlib_dd.hpp"
#include "gnstlib_polyeval.hpp"
#include "gnstlib_utils.hpp"

typedef std::complex<double> cmplx;

// Ei series around zero of Ei(x)
double ei_zero(const double x)
{
  const double x0 = 3.7250741078136663462e-1;
  const double a4 = 0.0021628544957325975355086703315891;
  const double a3 = 0.0132527785464987862931359100770888;
  const double a2 = 0.0687123596149106974230884712959543;
  const double a1 = 0.296094182264369238467577091129603;
  const double a0 = 1.10140810590373341810580422170947;

  double y = x - x0;
  dd s;

  // analytic continuation 2F2 with z0 = 0.3725...
  dd_mul_dd_d(y, a4, s);
  dd_add_dd_d_ip(s, a3);
  dd_mul_dd_d_ip(s, y);
  dd_add_dd_d_ip(s, a2);
  dd_mul_dd_d_ip(s, y);
  dd_add_dd_d_ip(s, a1);
  dd_mul_dd_d_ip(s, y);
  dd_add_dd_d_ip(s, a0);

  dd_mul_dd_d_ip(s, x);
  dd_add_dd_d_ip(s, GNSTLIB::constants::eulmasc + std::log(x));
  return s.hi + s.lo;
}

// Exponential integral Ei for real argument
double GNSTLIB::ei(const double x, int& err_id)
{
  double ans, frac, sump, sumq, t, xx0, xmx0, y, y2, w;
  double px[10], qx[10];

  const double xmax = 716.351;
  const double x0   = 3.7250741078136663466e-1;
  const double x01  = 381.5;
  const double x11  = 1024.0;
  const double x02  = -5.1182968633365538008e-5;

  // special values
  if (x == 0.0) {
    err_id = 4;
    return -inf;
  }
  else if (x == inf || x >= xmax) {
    err_id = 1;
    return inf;
  }
  else if (x == -inf) {
    err_id = 2;
    return 0.0;
  }

  // series around zero: 3.7250741078136663466e-1;
  if (x > 0.372505 && x < 0.372509)
    return ei_zero(x);

  if (x < 0.0) {
    y = std::abs(x);
    if (y < 1.0) {
      double P[7] = {
       -1.4815102102575750838086e5, 1.5026059476436982420737e5,
        8.9904972007457256553251e4, 1.5924175980637303639884e4,
        2.1500672908092918123209e3, 1.1669552669734461083368e2,
        5.0196785185439843791020e0
      };

      double Q[7] = {
        2.5666493484897117319268e5, 1.8434070063353677359298e5,
        5.2440529172056355429883e4, 8.1258035174768735759855e3,
        7.5043163907103936624165e2, 4.0205465640027706061433e1,
        1.000000000000000000000000,
      };

      ans = std::log(y) - polyeval_unroll7(P, y) / polyeval_unroll7(Q, y);
    } 
    else if (y <= 4.0) {

      double P[9] = {
        1.737331760720576030932e-8, 9.999989642347613068437e-1,
        1.487967702840464066613e+1, 7.633628843705946890896e+1,
        1.698106763764238382705e+2, 1.700632978311516129328e+2,
        7.246689782858597021199e+1, 1.107326627786831743809e+1,
        3.828573121022477169108e-1,
      };

      double Q[9] = {
        1.000000000000000000000e+0, 1.587964570758947927903e+1,
        9.021658450529372642314e+1, 2.342573504717625153053e+2,
        2.953136335677908517423e+2, 1.775728186717289799677e+2,
        4.662179610356861756812e+1, 4.344836335509282083360e+0,
        8.258160008564488034698e-2,
      };

      w = 1.0 / y;
      ans = -polyeval_unroll9(P, w) * std::exp(-y) / polyeval_unroll9(Q, w);
    } 
    else {

      double P[10] = {
        9.9999999999999999087819e-1, 5.2199632588522572481039e+1,
        1.0611777263550331766871e03, 1.0816852399095915622498e+4,
        5.9346841538837119172356e+4, 1.7503273087497081314708e+5,
        2.6181454937205639647381e+5, 1.7283375773777593926828e+5,
        3.5846198743996904308695e+4, 1.3276881505637444622987e+2,
      };

      double Q[10] = {
        1.000000000000000000000e+0,  5.4199632588522559414924E+1,
        1.1635769915320848035459E+3, 1.2842808586627297365998E+4,
        7.9231787945279043698718E+4, 2.7858134710520842139357E+5,
        5.4616842050691155735758E+5, 5.5903756210022864003380E+5,
        2.5989762083608489777411E+5, 3.9147856245556345627078E+4,
      };
      w = 1.0 / y;
      ans = -w * std::exp(-y) * (1.0 - w * 
        polyeval_unroll10(P, w) / polyeval_unroll10(Q, w));
    }
  } 
  else if (x < 6.0) {
    t = 0.6666666666666666666667 * x - 2.0;

    double P[10] = {
      -1.2963702602474830028590E01, -1.2831220659262000678155E03,
      -1.4287072500197005777376E04, -1.4299841572091610380064E06,
      -3.1398660864247265862050E05, -3.5377809694431133484800E08,
      3.1984354235237738511048E08,  -2.5301823984599019348858E10,
      1.2177698136199594677580E10,  -2.0829040666802497120940E11
    };

    double Q[10] = {
      7.6886718750000000000000E01, -5.5648470543369082846819E03,
      1.9418469440759880361415E05, -4.2648434812177161405483E06,
      6.4698830956576428587653E07, -7.0108568774215954065376E08,
      5.4229617984472955011862E09, -2.8986272696554495342658E10,
      9.8900934262481749439886E10, -8.9673749185755048616855E10
    };

    px[0] = 0.0;
    qx[0] = 0.0;
    px[1] = P[0];
    qx[1] = Q[0];
    for (int k = 1; k <= 8; k++) {
      px[k + 1] = t * px[k] - px[k - 1] + P[k];
      qx[k + 1] = t * qx[k] - qx[k - 1] + Q[k];
    }
    sump = 0.5 * t * px[9] - px[8] + P[9];
    sumq = 0.5 * t * qx[9] - qx[8] + Q[9];
    frac = sump / sumq;
    xmx0 = (x - x01 / x11) - x02;

    if (std::abs(xmx0) >= 0.037) {
      ans = std::log(x / x0) + xmx0 * frac;
    } 
    else {

      double P[4] = {
        3.5687548468071500413E+02, -5.4989956895857911039E+02,
        2.3642701335621505212E+02, -2.4562334077563243311E+01,
      };

      double Q[4] = {
        1.7843774234035750207E+02, -3.3442903192607538956E+02,
        1.9400230218539473193E+02, -3.5553900764052419184E+01,
      };

      xx0 = x + x0;
      y = xmx0 / (xx0);
      y2 = y * y;
      sump = polyeval_unroll4(P, y2);
      sumq = polyeval_unroll4(Q, y2);

      ans = (sump / (sumq * (xx0)) + frac) * xmx0;
    }
  } 
  else if (x < 12.0) {

    double R[10] = {
      -2.645677793077147237806E00, -2.378372882815725244124E00,
      -2.421106956980653511550E01,  1.052976392459015155422E01,
      1.945603779539281810439E01,  -3.015761863840593359165E01,
      1.120011024227297451523E01,  -3.988850730390541057912E00,
      9.565134591978630774217E00,   9.981193787537396413219E-1
    };

    double S[9] = {
      1.598517957704779356479E-4,  4.644185932583286942650E00,
      3.697412299772985940785E02, -8.791401054875438925029E00,
      7.608194509086645763123E02,  2.852397548119248700147E01,
      4.731097187816050252967E02, -2.369210235636181001661E02,
      1.249884822712447891440E00
    };

    frac = 0.0;
    for (int k = 0; k <= 8; k++) {
      frac = S[k] / (R[k] + x + frac);
    }
    ans = std::exp(x) * (R[9] + frac) / x;
  } 
  else if (x < 24.0) {

    double P[10] = {
      -1.647721172463463140042E00, -1.860092121726437582253E01,
      -1.000641913989284829961E01, -2.105740799548040450394E01,
      -9.134835699998742552432E-1, -3.323612579343962284333E01,
      2.495487730402059440626E01,   2.652575818452799819855E01,
      -1.845086232391278674524E00,  9.999933106160568739091E-1
    };

    double Q[9] = {
      9.792403599217290296840E01, 6.403800405352415551324E01,
      5.994932325667407355255E01, 2.538819315630708031713E02,
      4.429413178337928401161E01, 1.192832423968601006985E03,
      1.991004470817742470726E02,-1.093556195391091143924E01,
      1.001533852045342697818E00
    };

    frac = 0.0;
    for (int k = 0; k <= 8; k++) {
      frac = Q[k] / (P[k] + x + frac);
    }
    ans = std::exp(x) * (P[9] + frac) / x;
  } 
  else {

    double P[10] = {
      1.75338801265465972390E02,  -2.23127670777632409550E02,
      -1.81949664929868906455E01, -2.79798528624305389340E01,
      -7.63147701620253630855E00, -1.52856623636929636839E01,
      -7.06810977895029358836E00, -5.00006640413131002475E00,
      -3.00000000320981265753E00,  1.00000000000000485503E00
    };

    double Q[9] = {
      3.97845977167414720840E04,   3.97277109100414518365E00,
      1.37790390235747998793E02,   1.17179220502086455287E02,
      7.04831847180424675988E01,  -1.20187763547154743238E01,
      -7.99243595776339741065E00, -2.99999894040324959612E00,
      1.99999999999048104167E00
    };

    y = 1.0 / x;
    frac = 0.0;
    for (int k = 0; k <= 8; k++) {
      frac = Q[k] / (P[k] + x + frac);
    }
    frac += P[9];
    ans = y + y * y * frac;
    if (x <= xmax - 24.0) {
      ans *= std::exp(x);
    } else {
      // reformulated to avoid premature overflow
      ans = (ans * std::exp(x - 40.0)) * 2.3538526683701998541e17;
    }
  }

  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// Exponential integral E1 for real argument
double GNSTLIB::e1(const double x, int& err_id)
{
  double ans;

  // TODO: special cases
  // infinity + - and 0

  ans = -GNSTLIB::ei(-x, err_id);
  // warning: only the real part is returned
  if (err_id == 0 && x < 0.0) err_id = 100;
  return ans;
}

// Logarithmic integral for real positive argument
double GNSTLIB::li(const double x, int& err_id)
{
  err_id = 0;

  // special case x < 0
  if (x < 0.0) {
    err_id = 4;
    return nan; // TODO: review this
  }

  // special cases
  if (x == 0.0)
    return 0.0; // no err_id
  else if (std::abs(x - GNSTLIB::constants::soldner) < epsilon)
    return 0.0; // zero of the Li(x) for x \in (0, inf)
  else if (x == 1.0) {
    err_id = 5;
    return -inf;
  }

  return GNSTLIB::ei(std::log(x), err_id);
}

// numerical methods for E1(z), z complex argument
cmplx e1_series(const cmplx z, int& err_id)
{
  const int maxiter = 100;

  int k;
  double d;
  cmplx n, s, sp, u;

  err_id = 0;

  u = -GNSTLIB::constants::eulmasc - std::log(z);

  n = cmplx(1.0);
  d = 1.0;
  s = cmplx(0.0);
  sp = s;

  for (k = 1; k <= maxiter; k++) {
    n *= -z;
    d *= k;
    s += n / (k * d);
    if (convergence(s, sp))
      return u - s;
    else
      sp = s;
  }

  err_id = -1; // no convergence - use backup method
  return nan;
}

// e1 laguerre series
cmplx e1_laguerre(const cmplx z, int& err_id)
{
  const int maxiter = 100;
  const double tol = epsilon * 0.1;

  int d, k;
  cmplx lk, lk1, q, r, s;

  lk = cmplx(1.0);
  lk1 = z + cmplx(1.0);

  s = 1.0 / lk1;
  d = 1;

  for (k = 1; k <= maxiter; k++) {
    d += 1;
    q = (z + 2.0 * k + 1.0) / (k + 1.0) * lk1 - k / (k + 1.0) * lk;
    r = 1.0 / (static_cast<double>(d) * (q * lk1));
    s += r;
    lk = lk1;
    lk1 = q;
    if (std::abs(r) < tol) {
      return s * std::exp(-z);
    }
  }
  err_id = -1; // no convergence
  return nan;
}

// Exponential integral E1 for complex argument
cmplx GNSTLIB::e1(const cmplx z, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  const double mod  = std::abs(z);

  cmplx ans;
  err_id = 0;

  // special cases
  if (mod == 0.0) // pole
    return cmplx(inf); // argument is undefined, return complex infinity

  // special case E1(-x)
  if (im_z == 0.0 && re_z < 0.0)
    return cmplx(-GNSTLIB::ei(-re_z, err_id), -GNSTLIB::constants::pi);

  if (mod <= 1.5 || (mod <= 5.0 && re_z < 0.0))
    ans = e1_series(z, err_id);
  else if (re_z < -2.0 * std::abs(im_z)) {
    if (mod < 50.0) {
      if (std::abs(im_z) > 5.0)
        ans = e1_laguerre(z, err_id);
      else 
        ans = e1_series(z, err_id);
    } 
    else if (im_z > 0.0)
      ans = e1_laguerre(z, err_id) - cmplx(0.0, GNSTLIB::constants::pi);
    else
      ans = e1_laguerre(z, err_id) + cmplx(0.0, GNSTLIB::constants::pi);
  }
  else 
    ans = e1_laguerre(z, err_id);

  if (err_id == -1)
    return nan;

  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// Exponential integral Ei for complex argument
cmplx GNSTLIB::ei(const cmplx z, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  const double mod  = std::abs(z);

  cmplx ans; 
  err_id = 0;

  // complex special cases
  if (mod == 0.0) {
    err_id = 4;
    return cmplx(-inf);
  }
  else if (re_z == 0.0) {
    if (im_z == inf)
      return cmplx(0.0, GNSTLIB::constants::pi);
    else if (im_z == -inf)
      return cmplx(0.0, -GNSTLIB::constants::pi);
  } 
  else if (mod == inf){
    err_id = (re_z > 0.0) ? 1 : 2;
    return cmplx(re_z);
  }

  // real case
  if (im_z == 0.0)
    return GNSTLIB::ei(re_z, err_id);

  // main algorithm in terms of E1(z)
  if (im_z > 0.0)
    return GNSTLIB::e1(-z, err_id) + cmplx(0.0, GNSTLIB::constants::pi);
  else
    return -GNSTLIB::e1(-z, err_id) - cmplx(0.0, GNSTLIB::constants::pi);
}

// Logarithmic integral li for negative and complex argument
cmplx GNSTLIB::li(const cmplx z, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);

  // real case, x >= 0
  if (im_z == 0.0 && re_z >= 0.0)
    return GNSTLIB::li(re_z, err_id);
  else
    return GNSTLIB::ei(std::log(z), err_id);
}

// Numerical inverse of Exponential integral Ei(x)
double GNSTLIB::invei(const double x, int& err_id)
{
  const double tol = 10.0 * epsilon; // for large arguments
  const int maxiter = 10;

  int k;
  double an, nn, un, xn, xn1;

  err_id = 0;

  // argument lifting
  if (x >= maxdouble) {
    err_id = 1;
    return inf;
  }

  // starting point
  if (x > 1.0e6)
    xn = std::log(x * std::log(x * std::log(x)));
  else if (x > 10.0)
    xn = std::log(x * std::log(x));
  else if (x > 1.0)
    xn = std::log(x) + GNSTLIB::constants::eulmasc;
  else if (x <= 1.0 && x > -0.5)
    xn = GNSTLIB::constants::eulmasc;
  else if (x <= -0.5 && x > -40.0)
    xn = std::exp(x - GNSTLIB::constants::eulmasc);
  else {
    xn = std::exp(x - GNSTLIB::constants::eulmasc);
    err_id = gnstlib_fp_error_handler(xn); 
    return xn;
  }

  // Halley iteration (2 - 3) iterations for most cases.
  for (k = 0; k <= maxiter; k++) {
    an = GNSTLIB::ei(xn, err_id) - x;
    un = an * std::exp(-xn);
    nn = un * xn * 1.0 / (1.0 - 0.5 * un * (xn - 1.0));
    xn1 = xn;
    xn -= nn;
    if (std::abs(xn - xn1) / xn < tol) {
      err_id = gnstlib_fp_error_handler(xn);
      return xn;
    }
  }
  err_id = -1;
  return nan;  
}

// Numerical inverse of Logarithmic integral li(x)
 double GNSTLIB::invli(const double x, int& err_id) 
 {
  const double tol1 = 10.0 * epsilon;
  const double tol2 = 1000.0 * epsilon;
  const int maxiter = 10;

  int k;
  double an, nn, tol, xn, xn1;
  
  err_id = 0;

  // lifting for large parameters
  if (x >= maxdouble / 1000.0 || x == inf){
    err_id = 1;
    return inf;
  }

  tol = (x > 1e10) ? tol2 : tol1;    
    
  // starting value
  if (x > 3.5)
    xn = x * std::log(x);
  else if (x <= 3.5 && x > 0.75)
    xn = 1.0 + x;
  else if (x <= 0.75 && x > -0.5)
    xn = GNSTLIB::constants::soldner + GNSTLIB::constants::eulmasc * x;
  else if (x <= -0.5 && x > -36)
    xn = 1.0 + std::exp(x - GNSTLIB::constants::eulmasc);
  else 
    return 1.0;

  // Halley iteration (2-3 iterations for most cases)
  for (k = 0; k <= maxiter; k++) {
    an = GNSTLIB::li(xn, err_id) - x;
    nn = an * std::log(xn) * 1.0 / (1.0 + 0.5 * an / xn);
    xn1 = xn;
    xn = xn - nn;
    if (std::abs((xn - xn1) / xn) < tol){
      err_id = gnstlib_fp_error_handler(xn);
      return xn;
    }
  }
  err_id = -1;
  return nan;  
 }