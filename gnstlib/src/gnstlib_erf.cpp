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
#include <vector>

#include "faddeeva.hpp"
#include "gnstlib.hpp"
#include "gnstlib_polyeval.hpp"
#include "gnstlib_utils.hpp"

typedef std::complex<double> cmplx;

// Error function for real argument
double GNSTLIB::erf(const double x, int& err_id)
{
  double ans;

  err_id = 0;

  if (x == 0.0)
    return 0.0;

  // no special err_ids

  ans = std::erf(x);

  // floating-point error handler
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// Complementary error function
double GNSTLIB::erfc(const double x, int& err_id)
{
  double ans;

  err_id = 0;

  if (x == inf)
    return 0.0;

  // no special err_ids

  ans = std::erfc(x);

  // floating-point error handler
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// Scaled complementary error function for real argument
double GNSTLIB::erfcx(const double x, int& err_id)
{
  double ans;

  err_id = 0;

  if (x == inf)
    return 0.0;

  ans = Faddeeva::erfcx(x);

  err_id = gnstlib_fp_error_handler(ans); 
  return ans;   
}

// Imaginary error function for real argument
double GNSTLIB::erfi(const double x, int& err_id)
{
  double ans; 

  err_id = 0;

  if (x == 0.0)
    return 0.0;

  ans = Faddeeva::erfi(x);

  err_id = gnstlib_fp_error_handler(ans); 
  return ans; 
}

// Dawson integral for real argument
double GNSTLIB::dawson(const double x, int& err_id)
{
  double ans; 

  ans = Faddeeva::Dawson(x);

  err_id = gnstlib_fp_error_handler(ans); 
  return ans; 
}

// Error function for complex argument
cmplx GNSTLIB::erf(const cmplx z, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  cmplx ans;

  // special cases and error handling
  if (re_z == 0.0 || std::isnan(std::abs(re_z))) {
    if (std::abs(im_z) == inf) {
      // note: python 0 + i*inf => (nan+infj)
      err_id = 1;
      return z;
    }
  }

  ans = Faddeeva::erf(z);

  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// Complementary error function for complex argument
cmplx GNSTLIB::erfc(const cmplx z, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  cmplx ans;

  // special cases and error handling
  if (re_z == 0.0 || std::isnan(std::abs(re_z))) {
    if (std::abs(im_z) == inf) {
      // note: python 0 + i*inf => (nan+infj)
      err_id = 1;
      return -z;
    }
  }

  ans = Faddeeva::erfc(z);

  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// Scaled complementary error function for complex argument
cmplx GNSTLIB::erfcx(const cmplx z, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  cmplx ans;

  // special cases and error handling
  if (re_z == 0.0 || std::isnan(std::abs(re_z))) {
    if (std::abs(im_z) == inf) {
      // note: python 0 + i*inf => (nan+infj)
      err_id = 2;
      return 0;
    }
  }

  ans = Faddeeva::erfcx(z);

  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// Imaginary error function for complex argument
cmplx GNSTLIB::erfi(const cmplx z, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  cmplx ans;

  // special cases and error handling
  if (re_z == 0.0 || std::isnan(std::abs(re_z))) {
    if (std::abs(im_z) == inf) {
      // note: python 0 + i*inf => (nan+infj)
      err_id = 0;
      return cmplx(0.0, std::copysign(1.0, im_z));
    }
  }

  ans = Faddeeva::erfi(z);
  
  err_id = gnstlib_fp_error_handler(ans); 
  return ans; 
}

// Dawson function for complex argument
cmplx GNSTLIB::dawson(const cmplx z, int& err_id)
{
  cmplx ans;

  ans = Faddeeva::Dawson(z);

  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// Faddeeva / scaled complex error function for complex argument
cmplx GNSTLIB::faddeeva(const cmplx z, int& err_id)
{
  cmplx ans;

  ans = Faddeeva::w(z);

  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

template<typename T>
T fc_series(const T z, int& err_id)
{
  int k, maxiter = 50;
  double t, q, w;
  T s, sp, u, z4;

  err_id = 0;

  z4 = z * z;
  z4 *= z4;
  w = 2.46740110027233965470862275; // (pi/2)**2

  u = z;
  t = 1.0;
  q = 1.0;
  s = u;
  sp = s;

  for (k = 1; k <= maxiter; k++) {
    u *= -z4 * w;
    t *= (2*k - 1) * (2*k);
    q += 4;
    s += u / (t * q);
    if (convergence(s, sp))
      return s;
    else
      sp = s;
  }
  err_id = -1;
  return nan;
}

template<typename T>
T fs_series(const T z, int& err_id)
{
  int k, maxiter = 50;
  double t, q, w;
  T s, sp, u, z2, z4;

  err_id = 0;

  z2 = z * z;
  z4 = z2 * z2;
  w = 2.46740110027233965470862275; // (pi/2)**2

  u = z * z2 * GNSTLIB::constants::pihalf;
  t = 1.0;
  q = 3.0;
  s = u / q;
  sp = s;

  for (k = 1; k <= maxiter; k++) {
    u *= -z4 * w;
    t *= (2*k) * (2*k + 1);
    q += 4;
    s += u / (t * q);
    if (convergence(s, sp))
      return s;
    else
      sp = s;
  }
  err_id = -1;
  return nan;
}

void fg_rational(const double x, double &f, double &g)
{
  double u, t, x2;

  x2 = x * x;
  t = GNSTLIB::constants::pi * x2;
  u = 1.0 / (t * t);
  t = 1.0 / t;

  double PF[10] = {
    3.76329711269987889006e-20, 1.34283276233062758925e-16,
    1.72010743268161828879e-13, 1.02304514164907233465e-10,
    3.05568983790257605827e-8,  4.63613749287867322088e-6,
    3.45017939782574027900e-4,  1.15220955073585758835e-2,
    1.43407919780758885261e-1,  4.21543555043677546506e-1
  };

  double QF[11] = {
    1.25443237090011264384e-20, 4.52001434074129701496e-17,
    5.88754533621578410010e-14, 3.60140029589371370404e-11,
    1.12699224763999035261e-8,  1.84627567348930545870e-6,
    1.55934409164153020873e-4,  6.44051526508858611005e-3,
    1.16888925859191382142e-1,  7.51586398353378947175e-1,
    1.0
  };

  double PG[11] = {
    1.86958710162783235106e-22, 8.36354435630677421531e-19,
    1.37555460633261799868e-15, 1.08268041139020870318e-12,
    4.45344415861750144738e-10, 9.82852443688422223854e-8,
    1.15138826111884280931e-5,  6.84079380915393090172e-4,
    1.87648584092575249293e-2,  1.97102833525523411709e-1,
    5.04442073643383265887e-1
  };

  double QG[12] = {
    1.86958710162783236342E-22, 8.39158816283118707363E-19,
    1.38796531259578871258E-15, 1.10273215066240270757E-12,
    4.60680728146520428211E-10, 1.04314589657571990585E-7,
    1.27545075667729118702E-5,  8.14679107184306179049E-4,
    2.53603741420338795122E-2,  3.37748989120019970451E-1,
    1.47495759925128324529E0,   1.0
  };

  f = 1.0 - u * polyeval_unroll10(PF, u) / polyeval_unroll11(QF, u);
  g = t * polyeval_unroll11(PG, u) / polyeval_unroll12(QG, u);
}

void fg_asymptotic(const double x, double &f, double &g)
{
  int k, maxiter = 50;
  double d, fs, fsp, gs, gsp, q, r, w, t, u, y;

  r = 1.0 / (GNSTLIB::constants::pi * x);
  w = 1.0 / (GNSTLIB::constants::pihalf * x * x);
  q = w;
  w *= w;
  u = 1.0;
  d = 1.0;

  fs = 1.0;
  gs = 0.5;
  fsp = fs;
  gsp = gs;

  for (k = 1; k <= maxiter; k++) {
    y = 0.5 + 2*k;
    u *= -(y - 1) * (y - 2);
    t = u * y;
    d *= w;
    fs += u * d;
    gs += t * d;
    if (convergence(fs, fsp) && convergence(gs, gsp))
      break;
    else{
      fsp = fs; 
      gsp = gs;
    }
  }
  
  f = r * fs;
  g = r * q * gs;
}

double fc_rational(const double x)
{
  double ans, c, f, g, s, x2, x22, x4;

  x2 = x * x;
  if (x2 < 2.5625) {
    x4 = x2 * x2;
    double P[6] = {
      9.99999999999999998822e-1, -2.05525900955013891793e-1,
      1.88843319396703850064e-2, -6.45191435683965050962e-4, 
      9.50428062829859605134e-6, -4.98843114573573548651e-8
    };

    double Q[7] = {
      1.00000000000000000118e0,  4.12142090722199792936e-2,
      8.68029542941784300606e-4, 1.22262789024179030997e-5,
      1.25001862479598821474e-7, 9.15439215774657478799e-10,
      3.99982968972495980367e-12 
    };

    ans = x * polyeval_unroll6(P, x4) / polyeval_unroll7(Q, x4);
  }
  else {
    x22 = x2 * 0.5;
    c = GNSTLIB::cospi(x22);
    s = GNSTLIB::sinpi(x22);
    
    fg_rational(x, f, g);
    ans = 0.5 + (f * s - g * c) / (GNSTLIB::constants::pi * x);
  }
  return ans;
}

double fs_rational(const double x)
{
  double ans, c, f, g, s, x2, x22, x4;

  x2 = x * x;
  if (x2 < 2.5625) {
    x4 = x2 * x2;
    double P[6] = {
      3.18016297876567817986e11, -4.42979518059697779103e10,
      2.54890880573376359104e9,  -6.29741486205862506537e7,
      7.08840045257738576863e5,  -2.99181919401019853726e3
    };

    double Q[7] = {
      6.07366389490084639049e11, 2.24411795645340920940e10,
      4.19320245898111231129e8,  5.17343888770096400730e6,
      4.55847810806532581675e4,  2.81376268889994315696e2,
      1.0,
    };

    ans = x * x2 * polyeval_unroll6(P, x4) / polyeval_unroll7(Q, x4);
  }
  else {
    x22 = x2 * 0.5;
    c = GNSTLIB::cospi(x22);
    s = GNSTLIB::sinpi(x22);
    
    fg_rational(x, f, g);
    ans = 0.5 - (f * c + g * s) / (GNSTLIB::constants::pi * x);
  }

  return ans;
}

double GNSTLIB::fresnelc(const double x, int& err_id)
{
  double ans, c, f, g, s, w, z;

  z = std::abs(x);

  err_id = 0;

  if (z == 0.0)
    return 0.0;
  if (z == inf)
    return 0.5 * sign(x);
  
  if (z <= 1.0)
    ans = fc_series(z, err_id);
  else if (z < 5.5)
    ans = fc_rational(z);
  else {
    fg_asymptotic(z, f, g);
    w = z * z * 0.5;
    ans = 0.5 + f * GNSTLIB::sinpi(w) - g * GNSTLIB::cospi(w);
  }

  err_id = gnstlib_fp_error_handler(ans);
  return ans * sign(x);
}

double GNSTLIB::fresnels(const double x, int& err_id)
{
  double ans, c, f, g, s, w, z;

  z = std::abs(x);

  err_id = 0;

  if (z == 0.0)
    return 0.0;
  if (z == inf)
    return 0.5 * sign(x);
  
  if (z <= 1.0)
    ans = fs_series(z, err_id);
  else if (z < 5.5)
    ans = fs_rational(z);
  else {
    fg_asymptotic(z, f, g);
    w = z * z * 0.5;
    ans = 0.5 - f * GNSTLIB::cospi(w) - g * GNSTLIB::sinpi(w);
  }

  err_id = gnstlib_fp_error_handler(ans);
  return ans * sign(x);
}

// Fresnel integral C for complex argument
cmplx GNSTLIB::fresnelc(const cmplx z, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  cmplx erfi1, erfi2, t, u, v;

  err_id = 0;

  // special cases and error handling
  if (re_z == 0.0 || std::isnan(std::abs(re_z))) {
    if (std::abs(im_z) == inf) {
      // note: python 0 + i*inf => (nan+infj)
      return cmplx(0.0, 0.5 * sign(im_z));
    }
  }

  if (std::abs(z) < 1.0) {
    t = fc_series(z, err_id);
  }
  else {
    u = cmplx(1.0,-1.0);
    v = cmplx(1.0, 1.0);

    t = 0.5 * GNSTLIB::constants::sqrtpi * z;
    erfi1 = GNSTLIB::erfi(u * t, err_id);
    erfi2 = GNSTLIB::erfi(v * t, err_id);

    t = v * 0.25;
    t *= erfi1 - cmplx(0.0, 1.0) * erfi2; 
  }
  
  err_id = gnstlib_fp_error_handler(t);
  return t;
}

// Fresnel integral S for complex argument
cmplx GNSTLIB::fresnels(const cmplx z, int& err_id)
{
  const double re_z = std::real(z);
  const double im_z = std::imag(z);
  cmplx erfi1, erfi2, t, u, v;

  err_id = 0;

  // special cases and error handling
  if (re_z == 0.0 || std::isnan(std::abs(re_z))) {
    if (std::abs(im_z) == inf) {
      // note: python 0 + i*inf => (nan+infj)
      return cmplx(0.0, -0.5 * sign(im_z));
    }
  }

  if (std::abs(z) < 1.0) {
    t = fs_series(z, err_id);
  }
  else {
    u = cmplx(1.0,-1.0);
    v = cmplx(1.0, 1.0);

    t = 0.5 * GNSTLIB::constants::sqrtpi * z;
    erfi1 = GNSTLIB::erfi(u * t, err_id);
    erfi2 = GNSTLIB::erfi(v * t, err_id);

    t = cmplx(-1.0, 1.0) * 0.25;
    t *= erfi1 + cmplx(0.0, 1.0) * erfi2;
  }

  err_id = gnstlib_fp_error_handler(t);
  return t;
}

// Voigt profile
double GNSTLIB::voigt_profile(const double x, const double sigma, 
  const double gamma, int& err_id)
{
  double q, voigt;
  cmplx u, z;

  // check inputs
  if (sigma <= 0.0 || gamma <= 0.0) {
    err_id = 4;
    return nan;
  }

  // sigma = alpha / sqrt(2.0 * GNSTLIB::constants::log_2);
  q = sigma * GNSTLIB::constants::sqrttwopi;
  z = cmplx(x, gamma) / (sigma * GNSTLIB::constants::sqrt2);

  u = GNSTLIB::faddeeva(z, err_id);
  voigt = std::real(u) / q;

  err_id = gnstlib_fp_error_handler(voigt);
  return voigt;
}


// Inverse error function using rational approximants
// reference:
//  J. M. Blair, C. A. Edwards, and J. H. Johnson, "Rational Chebyshev 
//   approximations for the inverse of the error function," Math. Comp. 30,
//   pp. 827--830 (1976).
//  http://dx.doi.org/10.1090/S0025-5718-1976-0421040-7 
//  http://www.jstor.org/stable/2005402
double GNSTLIB::inverf(const double x, int& err_id)
{
  double a, ans, t;

  a = std::abs(x);

  // special cases
  if (a > 1.0) {
    err_id = 4;
    return nan;
  }

  if (x == 1.0) {
    err_id = 1;
    return inf;
  }

  if (x == -1.0) {
    err_id = 2;
    return -inf;
  }

  if (a <= 0.75) {
    t = x * x - 0.5625;

    double P[7] = {
      0.160304955844066229311e2, -0.90784959262960326650e2,
      0.18644914861620987391e3,  -0.16900142734642382420e3,
      0.6545466284794487048e2,   -0.864213011587247794e1,
      0.1760587821390590
    };

    double Q[7] = {
      0.147806470715138316110e2, -0.91374167024260313936e2,
      0.21015790486205317714e3,  -0.22210254121855132366e3,
      0.10760453916055123830e3,  -0.206010730328265443e2,
      0.1e1
    };

    ans = x * polyeval_unroll7(P, t) / polyeval_unroll7(Q, t);
  } 
  else if (a <= 0.9375) {
    t = x*x - 0.87890625;

    double P[8] = {
      -0.152389263440726128e-1,  0.3444556924136125216,
      -0.29344398672542478687e1, 0.11763505705217827302e2,
      -0.22655292823101104193e2, 0.19121334396580330163e2,
      -0.5478927619598318769e1,  0.237516689024448
    };

    double Q[8] = {
      -0.108465169602059954e-1,  0.2610628885843078511,
      -0.24068318104393757995e1, 0.10695129973387014469e2,
      -0.23716715521596581025e2, 0.24640158943917284883e2,
      -0.10014376349783070835e2, 0.1e1
    };

    ans = x * polyeval_unroll8(P, t) / polyeval_unroll8(Q, t);
  }
  else {
    t = 1.0 / std::sqrt(-std::log(1.0 - a));

    double P[9] = {
      0.10501311523733438116e-3,   0.1053261131423333816425e-1,
      0.26987802736243283544516,   0.23268695788919690806414e1,
      0.71678547949107996810001e1, 0.85475611822167827825185e1,
      0.68738088073543839802913e1, 0.3627002483095870893002e1,
      0.886062739296515468149
    };

    double Q[10] = {
      0.10501266687030337690e-3,    0.1053286230093332753111e-1,
      0.27019862373751554845553,    0.23501436397970253259123e1,
      0.76078028785801277064351e1,  0.111815861040569078273451e2,
      0.119487879184353966678438e2, 0.81922409747269907893913e1,
      0.4099387907636801536145e1,   0.1e1
    };

    ans = polyeval_unroll9(P, t) / 
      (std::copysign(t, x) * polyeval_unroll10(Q, t));
  }

  // floating-point error handler
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// inverse complementary error function using rational approximants
// based on julia implementation
double GNSTLIB::inverfc(const double x, int& err_id)
{
  double ans, t;

  if (x > 0.0625)
    return GNSTLIB::inverf(1.0 - x, err_id);

  // special cases
  if (x == 0.0) {
    err_id = 1;
    return inf;
  }

  if (x < 0.0) {
    err_id = 4;
    return nan;
  }

  if (x >= 1.0e-100) {
    t = 1.0 / std::sqrt(-std::log(x));

    double P[9] = {
      0.10501311523733438116e-3,   0.1053261131423333816425e-1,
      0.26987802736243283544516,   0.23268695788919690806414e1,
      0.71678547949107996810001e1, 0.85475611822167827825185e1,
      0.68738088073543839802913e1, 0.3627002483095870893002e1,
      0.886062739296515468149
    };

    double Q[10] = {
      0.10501266687030337690e-3,    0.1053286230093332753111e-1,
      0.27019862373751554845553,    0.23501436397970253259123e1,
      0.76078028785801277064351e1,  0.111815861040569078273451e2,
      0.119487879184353966678438e2, 0.81922409747269907893913e1,
      0.4099387907636801536145e1,   0.1e1
    };

    ans = polyeval_unroll9(P, t) / (t * polyeval_unroll10(Q, t));
  }
  else {
    t = 1.0 / std::sqrt(-std::log(x));

    double P[8] = {
      0.34654298588086350177e-9,    0.2508467920240757052055e-6,
      0.47378131963728602986534e-4, 0.313126037597786964083388e-2,
      0.77948764544143536994854e-1, 0.70045681233581643868271e0,
      0.18710420342167931668683e1,  0.714525477431351454283e0 
    };

    double Q[9] = {
      0.34654295673159511156e-9,    0.2508469079758802711487e-6,
      0.47379531295974913536339e-4, 0.313206353646177688480813e-2,
      0.78073489062764897214733e-1, 0.70715044799533758619993e0,
      0.19998515434911215105214e1,  0.1507290269273168000856e1,
      0.1e1
    };
  
    ans = polyeval_unroll8(P, t) / (t * polyeval_unroll9(Q, t));
  }

  // floating-point error handler
  err_id = gnstlib_fp_error_handler(ans); 
  return ans; 
}