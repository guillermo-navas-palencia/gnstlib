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
 * References:
 *  + A. J. MacLeod (1996b) Rational approximations, software and test methods 
 *    for sine and cosine integrals. Numer. Algorithms 12 (3-4), pp. 259–272. 
 */

#include <cmath>
#include <complex>
#include <vector>

#include "gnstlib.hpp"
#include "gnstlib_polyeval.hpp"
#include "gnstlib_utils.hpp"

typedef std::complex<double> cmplx;

// Sine integral Si for real argument
double GNSTLIB::si(const double x, int& err_id)
{
  double ans, f, g, z, z2;

  int si_sign = (x < 0.0) ? -1 : 1;

  err_id = 0;

  // special cases
  if (x == 0.0)
    return 0.0;
  else if (std::abs(x) == inf)
    return si_sign * GNSTLIB::constants::pihalf;

  z = std::abs(x);

  // other special cases
  if (z < 4.47e-8) {
    ans = z * si_sign;
    err_id = gnstlib_fp_error_handler(ans); 
    return ans;
  } else if (z > 2.32472e8) {
    if (z > 9.0072e15)
      ans = GNSTLIB::constants::pihalf;
    else {
      z2 = 1.0 / (z * z); 
      ans = GNSTLIB::constants::pihalf - std::cos(z) / z - std::sin(z) * z2;
    }
    ans *= si_sign;
    err_id = gnstlib_fp_error_handler(ans);
    return ans;
  }

  if (z <= 6.0) { 
    double P[8] = {
      1.00000000000000000000000,  -0.44663998931312457298e-1, 
      0.11209146443112369449e-2,  -0.13276124407928422367e-4, 
      0.85118014179823463879e-7,  -0.29989314303147656479e-9, 
      0.55401971660186204711e-12, -0.42406353433133212926e-15
    };

    double Q[8] = {
      1.00000000000000000000000,  0.10891556624243098264e-1, 
      0.59334456769186835896e-4,  0.21231112954641805908e-6,
      0.54747121846510390750e-9,  0.10378561511331814674e-11, 
      0.13754880327250272679e-14, 0.10223981202236205703e-17
    };

    z2 = z * z;
    ans = z * polyeval_unroll8(P, z2) / polyeval_unroll8(Q, z2);
  } 
  else if (z > 6.0 && z <= 12.0) {
    double PF[8] = {
      0.99999999962173909991,    0.36451060338631902917e3, 
      0.44218548041288440874e5,  0.22467569405961151887e7, 
      0.49315316723035561922e8,  0.43186795279670283193e9, 
      0.11847992519956804350e10, 0.45573267593795103181e9
    };

    double QF[8] = {
      1.0000000000000000000000,  0.36651060273229347594e3,
      0.44927569814970692777e5,  0.23285354882204041700e7,
      0.53117852017228262911e8,  0.50335310667241870372e9,
      0.16575285015623175410e10, 0.11746532837038341076e10
    };

    double PG[9] = {
      0.99999999920484901956,    0.51385504875307321394e3, 
      0.92293483452013810811e5,  0.74071341863359841727e7, 
      0.28142356162841356551e9,  0.49280890357734623984e10, 
      0.35524762685554302472e11, 0.79194271662085049376e11, 
      0.17942522624413898907e11
    };

    double QG[9] = {
      1.0000000000000000000000,  0.51985504708814870209e3,
      0.95292615508125947321e5,  0.79215459679762667578e7,
      0.31977567790733781460e9,  0.62273134702439012114e10,
      0.54570971054996441467e11, 0.18241750166645704670e12,
      0.15407148148861454434e12
    };

    z2 = 1.0 / (z * z);
    f = polyeval_unroll8(PF, z2) / (z * polyeval_unroll8(QF, z2));
    g = z2 * polyeval_unroll9(PG, z2) / polyeval_unroll9(QG, z2);
    ans = GNSTLIB::constants::pihalf - f * std::cos(z) - g * std::sin(z);
  }
  else {
    double PF[8] = {
      0.19999999999999978257e1, 0.22206119380434958727e4, 
      0.84749007623988236808e6, 0.13959267954823943232e9,
      0.10197205463267975592e11, 0.30229865264524075951e12, 
      0.27504053804288471142e13, 0.21818989704686874983e13
    };

    double QF[8] = {
      1.0000000000000000000000,  0.11223059690217167788e4,
      0.43685270974851313242e6,  0.74654702140658116258e8,
      0.58580034751805687471e10, 0.20157980379272098841e12,
      0.26229141857684496445e13, 0.87852907334918467516e13
    };

    double PG[9] = {
      0.59999999999999993089e1,  0.96527746044997139158e4, 
      0.56077626996568834185e7,  0.15022667718927317198e10, 
      0.19644271064733088465e12, 0.12191368281163225043e14, 
      0.31924389898645609533e15, 0.25876053010027485934e16, 
      0.12754978896268878403e16
    };

    double QG[9] = {
      1.0000000000000000000000,  0.16287957674166143196e4,
      0.96636303195787870963e6,  0.26839734750950667021e9,
      0.37388510548029219241e11, 0.26028585666152144496e13,
      0.85134283716950697226e14, 0.11304079361627952930e16,
      0.42519841479489798424e16
    };

    z2 = 1.0 / (z * z);
    f = (1.0 - z2 * polyeval_unroll8(PF, z2) / polyeval_unroll8(QF, z2)) / z;
    g = (1.0 - z2 * polyeval_unroll9(PG, z2) / polyeval_unroll9(QG, z2)) * z2;
    ans = GNSTLIB::constants::pihalf - f * std::cos(z) - g * std::sin(z);
  }

  ans *= si_sign;
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;  
}

template<typename T>
T si_series(const T z, int& err_id)
{
  int k, maxiter = 50;
  double d, m, r;
  T n, s, sp, t, u;

  err_id = 0;

  u = 0.25 * z * z;
  n = 1.0;
  r = 1.0;
  d = 1.0;
  m = 1.0;
  s = 1.0;
  sp = s;

  for (k = 1; k <= maxiter; k++) {
    n *= -u;
    d *= (0.5 + static_cast<double>(k));
    r *= (static_cast<double>(k) - 0.5);
    m *= k;
    t = n * r / (m * d * d);
    s += t;
    if (convergence(s, sp))
      return z * s;
    else
      sp = s;
  }
  err_id = -1;
  return nan;
}

// Hyperbolic sine integral Shi - series expansion - real and complex argument
template<typename T>
T shi_series(const T z, int& err_id)
{
  int k, maxiter = 50;
  double d, m, r;
  T n, s, sp, t, u;

  err_id = 0;

  u = 0.25 * z * z;
  n = 1.0;
  r = 1.0;
  d = 1.0;
  m = 1.0;
  s = 1.0;
  sp = s;

  for (k = 1; k <= maxiter; k++) {
    n *= u;
    d *= (0.5 + static_cast<double>(k));
    r *= (static_cast<double>(k) - 0.5);
    m *= k;
    t = n * r / (m * d * d);
    s += t;
    if (convergence(s, sp))
      return z * s;
    else
      sp = s;
  }
  err_id = -1;
  return nan;
}

// Hyperbolic sine integral Shi for real argument
double GNSTLIB::shi(const double x, int& err_id)
{
  double ans, z;

  int shi_sign = (x < 0.0) ? -1 : 1;

  // special cases
  if (x == 0.0)
    return 0.0;
  else if (std::abs(x) == inf) {
    err_id = 1;
    return x;
  }

  z = std::abs(x);
  
  if (z < 1.0)
    ans = shi_series(z, err_id);
  else {
    if (z < 20.0)
      ans = 0.5 * (GNSTLIB::ei(z, err_id) + GNSTLIB::e1(z, err_id));  
    else
      ans = 0.5 * GNSTLIB::ei(z, err_id); // e1(z) can be safely neglected.
  }

  err_id = gnstlib_fp_error_handler(ans);
  return ans;
}

// Cosine integral Ci for real argument
double GNSTLIB::ci(const double x, int& err_id)
{
  double ans, cx, cx2, dif, f, g, logval, sd, sn, sum, sx, z, z2;

  err_id = 0;

  z = std::abs(x);

  // special cases
  if (x == 0.0) {
    err_id = 4;
    return -inf;
  }

  // other special cases
  if (z <= 1.48996E-8)
    return GNSTLIB::constants::eulmasc + std::log(z);
  else if (z > 2.324953e8) {
    if (z > 1.4148475e16)
      ans = 0.0;  
    else {
      z2 = 1.0 / (z * z);
      ans = std::sin(z) / z - std::cos(z) * z2; 
    }
    err_id = gnstlib_fp_error_handler(ans);
    return ans;
  }

  if (z <= 3.0) {

    double P[6] = {
      -0.24607411378767540707,    0.72113492241301534559e-2, 
      -0.11867127836204767056e-3, 0.90542655466969866243e-6, 
      -0.34322242412444409037e-8, 0.51950683460656886834e-11
    };

    double Q[6] = {
      1.00000000000000000000000,  0.12670095552700637845e-1,
      0.78168450570724148921e-4,  0.29959200177005821677e-6,
      0.73191677761328838216e-9,  0.94351174530907529061e-12
    };

    z2 = z*z;
    sn = polyeval_unroll6(P, z2);
    sd = polyeval_unroll6(Q, z2);
    dif = (z - 0.6162109375) - 0.29454812071623379711e-3;
    sum = 0.6162109375 + 0.29454812071623379711e-3;
    if (std::abs(dif) < 0.046875) {
      cx = dif / (sum + z);
      cx2 = cx * cx;
      sx = 0.83930008362695945726e1 + cx2 * (-0.65306663899493304675e1 
        + cx2 * 0.569155722227490223);
      sx = sx / (0.41965004181347972847e1 + cx2 * (-0.46641666676862479585e1 
        + cx2));
      logval = cx * sx;
    } else {
      logval = std::log(z / sum);
    }
    ans = logval + dif * (z + sum) * sn / sd;
  }
  else if (z > 3.0 && z <= 6.0) {

    double P[8] = {
      -0.15684781827145408780,     0.66253165609605468916e-2,  
      -0.12822297297864512864e-3,  0.12360964097729408891e-5,
      -0.66450975112876224532e-8,  0.20326936466803159446e-10,
      -0.33590883135343844613e-13, 0.23686934961435015119e-16
    };

    double Q[7] = {
      1.0000000000000000000000,  0.96166044388828741188E-2,
      0.45257514591257035006E-4, 0.13544922659627723233E-6,
      0.27715365686570002081E-9, 0.37718676301688932926E-12,
      0.27706844497155995398E-15
    };

    z2 = z * z;
    sn = polyeval_unroll8(P, z2);
    sd = polyeval_unroll7(Q, z2);
    dif = (z - 3.3837890625) - 0.39136005118642639785e-3;
    sum = 3.3837890625 + 0.39136005118642639785e-3;
    if (std::abs(dif) < 0.25) {
      cx = dif / (sum + z);
      cx2 = cx * cx;
      sx = 0.83930008362695945726e1 + cx2 * (-0.65306663899493304675e1 
        + cx2 * 0.569155722227490223);
      sx = sx / (0.41965004181347972847e1 + cx2 * (-0.46641666676862479585e1 
        + cx2));
      logval = cx * sx;
    } else {
      logval = std::log(z / sum);
    }
    ans = logval + dif * (z + sum) * sn / sd;
  }
  else if (z > 6.0 && z <= 12.0) {

    double PF[8] = {
      0.99999999962173909991,    0.36451060338631902917e3, 
      0.44218548041288440874e5,  0.22467569405961151887e7, 
      0.49315316723035561922e8,  0.43186795279670283193e9, 
      0.11847992519956804350e10, 0.45573267593795103181e9
    };

    double QF[8] = {
      1.0000000000000000000000,  0.36651060273229347594e3,
      0.44927569814970692777e5,  0.23285354882204041700e7,
      0.53117852017228262911e8,  0.50335310667241870372e9,
      0.16575285015623175410e10, 0.11746532837038341076e10
    };

    double PG[9] = {
      0.99999999920484901956,    0.51385504875307321394e3, 
      0.92293483452013810811e5,  0.74071341863359841727e7, 
      0.28142356162841356551e9,  0.49280890357734623984e10, 
      0.35524762685554302472e11, 0.79194271662085049376e11, 
      0.17942522624413898907e11
    };

    double QG[9] = {
      1.0000000000000000000000,  0.51985504708814870209e3,
      0.95292615508125947321e5,  0.79215459679762667578e7,
      0.31977567790733781460e9,  0.62273134702439012114e10,
      0.54570971054996441467e11, 0.18241750166645704670e12,
      0.15407148148861454434e12
    };

    z2 = 1.0 / (z * z);
    f = polyeval_unroll8(PF, z2) / (z * polyeval_unroll8(QF, z2));
    g = z2 * polyeval_unroll9(PG, z2) / polyeval_unroll9(QG, z2);
    ans = f * std::sin(z) - g * std::cos(z);
  }
  else {
    double PF[8] = {
      0.19999999999999978257e1, 0.22206119380434958727e4, 
      0.84749007623988236808e6, 0.13959267954823943232e9,
      0.10197205463267975592e11, 0.30229865264524075951e12, 
      0.27504053804288471142e13, 0.21818989704686874983e13
    };

    double QF[8] = {
      1.0000000000000000000000,  0.11223059690217167788e4,
      0.43685270974851313242e6,  0.74654702140658116258e8,
      0.58580034751805687471e10, 0.20157980379272098841e12,
      0.26229141857684496445e13, 0.87852907334918467516e13
    };

    double PG[9] = {
      0.59999999999999993089e1,  0.96527746044997139158e4, 
      0.56077626996568834185e7,  0.15022667718927317198e10, 
      0.19644271064733088465e12, 0.12191368281163225043e14, 
      0.31924389898645609533e15, 0.25876053010027485934e16, 
      0.12754978896268878403e16
    };

    double QG[9] = {
      1.0000000000000000000000,  0.16287957674166143196e4,
      0.96636303195787870963e6,  0.26839734750950667021e9,
      0.37388510548029219241e11, 0.26028585666152144496e13,
      0.85134283716950697226e14, 0.11304079361627952930e16,
      0.42519841479489798424e16
    };
    z2 = 1.0 / (z * z);
    f = (1.0 - z2 * polyeval_unroll8(PF, z2) / polyeval_unroll8(QF, z2)) / z;
    g = (1.0 - z2 * polyeval_unroll9(PG, z2) / polyeval_unroll9(QG, z2)) * z2;
    ans = f * std::sin(z) - g * std::cos(z);
  }

  err_id = gnstlib_fp_error_handler(ans);

  // throw warning: only the real part is returned
  if (err_id == 0 && x < 0.0) err_id = 100;
  return ans;
}

template<typename T>
T ci_series(const T z, int& err_id) 
{
  int k, maxiter = 50;
  double t, q;
  T f, u, s, sp, z2;

  err_id = 0;

  z2 = z * z;
  u = 1.0;
  s = 0.0;
  t = 1.0;
  q = 0.0;
  sp = s;
  
  f = GNSTLIB::constants::eulmasc + std::log(z);

  for (k = 1; k <= maxiter; k++) {
    u *= -z2;
    t *= (2*k - 1) * (2*k);
    q += 2;
    s += u / (t * q);
    if (convergence(s, sp))
      return s + f;
    else 
      sp = s;
  }
  err_id = -1;
  return nan;
}

template<typename T>
T chi_series(const T z, int& err_id)
{
  int k, maxiter = 50;
  double d;
  T s, sp, t, u;

  err_id = 0;

  u = z * z;
  t = u;
  d = 2.0;
  s = t * 0.5;
  sp = s;

  for (k = 2; k <= maxiter; k++) {
    t *= u;
    d *= (2*k - 1) * (2*k);
    s += t / (k * d);
    if (convergence(s, sp))
      return std::log(z) + GNSTLIB::constants::eulmasc + 0.5 * s;
    else
      sp = s;
  }
  err_id = -1;
  return nan;
}

// Hyperbolic cosine Chi integral for real argument
double GNSTLIB::chi(const double x, int& err_id)
{
  double ans, z;

  // special cases
  if (x == 0.0) {
    err_id = 4;
    return -inf;
  } else if(std::abs(x) == inf) {
    err_id = 1;
    return x;
  }

  z = std::abs(x);

  if (z < 1.0) {
    if (z < 1e-20)
      ans = std::log(z) + GNSTLIB::constants::eulmasc;
    else
      ans = chi_series(z, err_id);
  }
  else {
    if (z < 20)
      ans = 0.5 * (GNSTLIB::ei(z, err_id) - GNSTLIB::e1(z, err_id));
    else
      ans = 0.5 * GNSTLIB::ei(z, err_id); // e1(z) can be safely neglected.
  }

  err_id = gnstlib_fp_error_handler(ans);

  if (err_id == 0 && x < 0.0) err_id = 100;
  return ans;
}

// Sine integral Si for complex argument
cmplx GNSTLIB::si(const cmplx z, int& err_id) 
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  const double mod  = std::abs(z);

  cmplx a, ans, y;
  cmplx i = cmplx(0.0, 1.0);

  // complex special cases
  if (mod == inf && re_z == 0.0) {
    err_id = 1;
    return z;
  }

  // real case
  if (im_z == 0.0)
    return GNSTLIB::si(re_z, err_id);

  if (mod < 4.0)
    return si_series(z, err_id);

  // Si(z) = 1/2*i*(E1(-iz) - E1(iz)) + pi/2, |ph z| < pi/2
  y = i * z;
  a = GNSTLIB::ei(y, err_id) - GNSTLIB::ei(-y, err_id);
  ans = -0.5 * i * a;
  if(re_z > 0.0)
    ans -= 0.5 * GNSTLIB::constants::pi;
  else if(re_z < 0.0)
    ans += 0.5 * GNSTLIB::constants::pi;

  err_id = gnstlib_fp_error_handler(ans);
  return ans;
}

// Hyperbolic sine integral Shi for complex argument
cmplx GNSTLIB::shi(const cmplx z, int& err_id) 
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  const double mod  = std::abs(z);
  
  cmplx ans;
  cmplx i = cmplx(0.0, 1.0);

  err_id = 0;

  // complex special cases
  if (re_z == 0.0 || std::isnan(std::abs(re_z))) {
    // note: python 0 + i*inf => (nan+infj)
    if (im_z == inf)
      return cmplx(0.0, GNSTLIB::constants::pihalf);
    else if (im_z == -inf)
      return cmplx(0.0, -GNSTLIB::constants::pihalf);
  }

  // real case
  if (im_z == 0.0)
    return GNSTLIB::shi(re_z, err_id);

  // Shi(z) = -i * Si(i * z)
  if (mod < 4.0)
    ans = shi_series(z, err_id);
  else
    ans = -i * GNSTLIB::si(i * z, err_id);

  err_id = gnstlib_fp_error_handler(ans);
  return ans;
}

// Cosine integral Ci for complex argument
cmplx GNSTLIB::ci(const cmplx z, int& err_id) 
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  const double mod  = std::abs(z);

  cmplx ans, y;
  cmplx i = cmplx(0.0, 1.0);

  // complex special cases
  if (std::abs(im_z) == inf && (re_z == 0.0 || std::isnan(std::abs(re_z)))) {
    // note: python 0 + i*inf => (nan+infj)
    err_id = 1;
    return inf;
  }

  // real case 
  if (im_z == 0.0) {
    if (re_z >= 0.0)
        return GNSTLIB::ci(re_z, err_id);
      else
        return cmplx(GNSTLIB::ci(mod, err_id), GNSTLIB::constants::pi);
  }
  
  if (mod < 4.0)
    return ci_series(z, err_id);

  y = i * z;
  ans = 0.5 * (GNSTLIB::ei(y, err_id) + GNSTLIB::ei(-y, err_id));
  if (re_z == 0.0) {
    if (im_z > 0.0)
      ans += i * GNSTLIB::constants::pihalf;
    else // z.imag < 0.0
      ans -= i * GNSTLIB::constants::pihalf;
  } else if (re_z < 0.0) {
    if (im_z > 0.0)
      ans += i * GNSTLIB::constants::pi;
    else
      ans -= i * GNSTLIB::constants::pi;
  }
  
  err_id = gnstlib_fp_error_handler(ans);
  return ans;
}

// Hyperbolic cosine integral Chi for complex argument
cmplx GNSTLIB::chi(const cmplx z, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  const double mod  = std::abs(z);
  
  cmplx ans;
  cmplx i = cmplx(0.0, 1.0);

  err_id = 0;

  // complex special cases
  if (re_z == 0.0 || std::isnan(std::abs(re_z))) {
    // note: python 0 + i*inf => (nan+infj)
    if (im_z == inf)
      return cmplx(0.0, GNSTLIB::constants::pihalf);
    else if (im_z == -inf)
      return cmplx(0.0, -GNSTLIB::constants::pihalf);
  }

  // real case
  if (im_z == 0.0) {
    if (re_z >= 0.0)
      return GNSTLIB::chi(re_z, err_id);
    else
      return cmplx(GNSTLIB::chi(-re_z, err_id), GNSTLIB::constants::pi);
  }

  // Chi(z) = Ci(i * z) + log(z) - log(i * z)
  if (std::abs(z) < 4.0)
    ans = chi_series(z, err_id);
  else
    ans = GNSTLIB::ci(i * z, err_id) + std::log(z) - std::log(i * z);

  err_id = gnstlib_fp_error_handler(ans);
  return ans;
}

// Compute sine and cosine integral at once for complex argument
void GNSTLIB::sici(const cmplx z, cmplx& si, cmplx& ci, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);

  cmplx _si, _ci, a, ei_p, ei_n, y;
  cmplx i = cmplx(0.0, 1.0);

  // TODO:
  // different flags for errors occurred in Si or Ci.

  // real case
  if (im_z == 0.0) {
    si = cmplx(GNSTLIB::si(re_z, err_id));
    if (re_z >= 0.0)
      ci = cmplx(GNSTLIB::ci(re_z, err_id));
    else
      ci = cmplx(GNSTLIB::ci(std::abs(re_z), err_id), GNSTLIB::constants::pi);
    return;
  }

  // special case |z| < 4 : call direct power series - faster
  if (std::abs(z) < 4.0) {
    si = si_series(z, err_id);
    ci = ci_series(z, err_id);
    return;
  }

  y = i * z;
  ei_p = GNSTLIB::ei(y, err_id);
  ei_n = GNSTLIB::ei(-y, err_id);

  _si = -0.5 * i * (ei_p - ei_n);
  _ci = 0.5 * (ei_p + ei_n);

  if (re_z == 0.0) {
    if (im_z > 0.0)
      _ci += i * GNSTLIB::constants::pihalf;
    else
      _ci -= i * GNSTLIB::constants::pihalf;
  }

  if (re_z > 0.0) {
    _si -= GNSTLIB::constants::pihalf;
  } else if (re_z < 0.0) {
    _si += GNSTLIB::constants::pihalf;
    if (im_z > 0.0)
      _ci += i * GNSTLIB::constants::pi;
    else
      _ci -= i * GNSTLIB::constants::pi;
  }
  si = _si;
  ci = _ci;
}