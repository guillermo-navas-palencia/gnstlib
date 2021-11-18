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
#include <quadmath.h>

#include "gnstlib.hpp"
#include "gnstlib_utils.hpp"

// dominant part is approximately of x^a * exp(-x) / gamma(a+1)
double dompart(const double a, const double x, const bool qt)
{
  double ans, c, dp, la, lambda, lnx, mu, r;
  int err_id;

  lnx = std::log(x);

  if (a <= 1.0)
    r = -x + a * lnx;
  else {
    if (x == a)
      r = 0.0;
    else {
      la = x / a;
      r = a * (1.0 - la + std::log(la));
    }
    r -= 0.5 * std::log(6.2832 * a);
  }

  if (r < explow)
    dp = 0.0;
  else
    dp = std::exp(r);

  if (qt)
    ans = dp;
  else {
    if (a < 8.0)
      ans = std::exp(a * lnx - x) / GNSTLIB::gamma(a + 1.0, err_id);
    else {
      ans = 1.0 / (std::sqrt(a * GNSTLIB::constants::twopi)
        * GNSTLIB::gammastar(a, err_id));
      lambda = x / a;
      if (lambda > 0.3 && lambda < 2.36) {
        mu = lambda - 1.0;
        c = GNSTLIB::log1pmx(mu);
        ans *= std::exp(a * c);
      } 
      else {
        ans *= std::exp(a * std::log(lambda) + a - x);
      }   
    }
  }
  return ans;
}

double alpha(const double x)
{
  if (x > 0.25)
    return x + 0.25;
  else if (x >= dwarf)
    return -0.6931 / std::log(x);
  else
    return -0.6931 / std::log(dwarf);
}

double ptaylor(const double a, const double x, const double dp, int& err_id)
{
  const int maxiter = 500;

  int k;
  double c, r, s, sp;

  err_id = 0;

  if (dp == 0.0)
    return 0.0;

  s = 1.0;
  sp = s;
  c = 1.0;
  r = a;
  for (k = 1; k <= maxiter; k++) {
    r += 1;
    c *= x / r;
    s += c;
    if (convergence(s, sp))
      return dp * s;
    else
      sp = s;
  }
  
  err_id = -1; // no convergence - use backup method
  return nan;
}

double qtaylor(const double a, const double x, const double dp, int& err_id)
{
  const int maxiter = 500;

  int k;
  double f, lnx, p, q, r, s, t, u, v, vp;

  err_id = 0;

  if (dp == 0.0)
    return 0.0;

  lnx = std::log(x);
  r = a * lnx;
  q = r * GNSTLIB::expm1x(r);
  s = a * (1.0 - a) * GNSTLIB::auxgam(a, err_id);
  q *= (1.0 - s);

  f = a * (1.0 - s) * std::exp((a + 1.0) * lnx) / (a + 1.0);

  u = s - q;
  p = a * x;
  q = a + 1.0;
  r = a + 3.0;
  t = 1.0;
  v = 1.0;
  vp = v;
  for (k = 1; k <= maxiter; k++) {
    p += x;
    q += r;
    r += 2;
    t *= -p/q;
    v += t;
    if (convergence(v, vp))
      return v * f + u;
    else
      vp = v;
  } 
  
  err_id = -1; // no convergence - use backup method
  return nan;
}

double qfraction(const double a, const double x, const double dp, int& err_id)
{
  const int maxiter = 500;

  int k;
  double f, g, gp, p, q, r, ro, s, t, tau, xp1ma;

  err_id = 0;

  if (dp == 0.0)
    return 0.0;

  s = 1.0 - a;
  xp1ma = x + s;
  q = (x - 1.0 - a) * xp1ma;
  r = 4.0 * xp1ma;
  f = (a / xp1ma) * dp;

  p = 0.0;
  ro = 0.0;
  t = 1.0;
  g = 1.0;
  gp = g;

  for (k = 1; k <= maxiter; k++)
  {
    p += s;
    q += r;
    r += 8.0;
    s += 2.0;
    tau = p * (1.0 + ro);
    ro = tau / (q - tau);
    t *= ro;
    g += t;
    if (convergence(g, gp))
      return g * f;
    else
      gp = g;
  }

  err_id = -1; // no convergence - use backup method
  return nan;
}

double saeta(const double a, const double eta)
{
  double s, t, y;
  std::vector<double> bm (27);
  std::vector<double> fm (27);

  fm[0] = 1.0;
  fm[1] = -1.0 / 3.0;
  fm[2] =  1.0 / 12.0;
  fm[3] = -2.0 / 135.0;
  fm[4] =  1.0 / 864.0;
  fm[5] = 1.0 / 2835.0;
  fm[6] = -139.0 / 777600.0;
  fm[7] =  1.0 / 25515.0;
  fm[8] = -571.0 / 261273600.0;
  fm[9] = -281.0 / 151559100.0;
  fm[10] =  8.29671134095308601e-7;
  fm[11] = -1.76659527368260793e-7;
  fm[12] =  6.70785354340149857e-9;
  fm[13] =  1.02618097842403080e-8;
  fm[14] = -4.38203601845335319e-9;
  fm[15] =  9.14769958223679023e-10;
  fm[16] = -2.55141939949462497e-11;
  fm[17] = -5.83077213255042507e-11;
  fm[18] =  2.43619480206674162e-11;
  fm[19] = -5.02766928011417559e-12;
  fm[20] =  1.10043920319561347e-13;
  fm[21] =  3.37176326240098538e-13;
  fm[22] = -1.39238872241816207e-13;
  fm[23] =  2.85348938070474432e-14;
  fm[24] = -5.13911183424257258e-16;
  fm[25] = -1.97522882943494428e-15;
  fm[26] =  8.09952115670456133e-16;

  bm[25] = fm[26];
  bm[24] = fm[25];
  for (int m = 24; m --> 1; ) {
    bm[m-1] = fm[m] + (m + 1) * bm[m + 1] / a;
  }
  s = bm[0];
  t = s;
  y = eta;

  for (int m = 1; m < 25; m++) {
    t = bm[m] * y;
    s += t;
    y *= eta;
    if (std::abs(t / s) < epsilon) {
      break;
    }
  }
  return s / (1.0 + bm[1] / a);
} 

double pqasymp(const double a, const double x, const double dp, const bool p)
{
  int err_id;
  double eta, mu, s, u, v, y;

  if (dp == 0.0) {
    if (p)
      return 0.0;
    else
      return 1.0;
  }

  s = (p) ? -1.0 : 1.0;

  mu = (x - a) / a;
  y = -GNSTLIB::log1pmx(mu);

  eta = (y < 0.0) ? 0.0 : std::sqrt(2.0 * y);
  y *= a;
  v = std::sqrt(std::abs(y));
  if (mu < 0.0) {
    eta = -eta;
    v = -v;
  }

  u = 0.5 * GNSTLIB::erfc(s * v, err_id);
  v = s * std::exp(-y) * saeta(a, eta) / sqrt(GNSTLIB::constants::twopi * a);
  return u + v;
}

void gammainc_ratios(const double a, const double x, double& p, double& q, 
  int& err_id)
{
  double dp, lnx;

  err_id = 0;

  if (a <= 0.0) {
    err_id = 4;
    p = nan; q = nan;
    return;
  }

  if (x < 0.0) {
    err_id = 5;
    p = nan; q = nan;
    return;
  }

  if (x < dwarf)
    lnx = std::log(dwarf);
  else
    lnx = std::log(x);

  if (a > alpha(x)) {
    dp = dompart(a, x, false);
    if (dp < 0.0) {
      err_id = 6;
      p = nan; q = nan;
      return;
    } 
    else {
      if (x < 0.34 * a || a < 16.5)
        p = ptaylor(a, x, dp, err_id); 
      else
        p = pqasymp(a, x, dp, true);
      q = 1.0 - p;
    }
  }
  else {
    if (a < -dwarf / lnx)
      q = 0.0;
    else {
      if (x < 1.0) {
        dp = dompart(a, x, true);
        if (dp < 0.0) {
          err_id = 6;
          p = nan; q = nan;
          return;
        }
        else {
          q = qtaylor(a, x, dp, err_id);
          p = 1.0 - q;
        }
      }
      else {
        dp = dompart(a, x, false);
        if (dp < 0.0) {
          err_id = 6;
          p = nan; q = nan;
          return;
        }
        else {
          if (x > 1.5 * a || a < 12.0)
            q = qfraction(a, x, dp, err_id);
          else {
            q = pqasymp(a, x, dp, false);
            if (dp == 0.0)
              q = 0.0;
          }
          p = 1.0 - q;
        }
      }
    }
  }
}

// Regularized lower incomplete gamma function P(a, x)
double GNSTLIB::gammainc_p(const double a, const double x, int& err_id)
{
  double p, q;

  err_id = 0;

  // special values
  if (x == 0.0) 
    return 0.0;

  gammainc_ratios(a, x, p, q, err_id);

  if (err_id == 0)
    err_id = gnstlib_fp_error_handler(p); 
    
  return p;
}

// Regularized upper incomplete gamma function Q(a, x)
double GNSTLIB::gammainc_q(const double a, const double x, int& err_id)
{
  double p, q;

  err_id = 0;

  // special cases 
  if (x == 0.0)
    return 1.0;
  else if (a == 1 && x > 0.0)
    return GNSTLIB::exp(-x);

  gammainc_ratios(a, x, p, q, err_id);

  if (err_id == 0)
    err_id = gnstlib_fp_error_handler(q);
    
  return q;
}

int use_expint_asymp_x(const double v, const double x)
{
  int i, n_max;
  double r, v1;

  n_max = static_cast<int>(ceil(x - v));
  v1 = v - 1.0;
  r = 1.0;
  for (i = 1; i <= n_max; i++) {
    r *= (v1 + i) / x;
    if (r < epsilon)
      return i;
  }
  return 0;
}

double expint_asymp_x(const double v, const double x, const int n) 
{
  int i;
  double d, s, u, v1;

  v1 = v - 1.0;

  // series
  u = 1.0;
  d = x;
  s = 1.0 / x;

  for (i = 1; i <= n; i++) {
    u *= -(v1 + i);
    d *= x;
    s += u / d;
  }
  return std::exp(-x) * s;
}

int use_expint_asymp_v(const double v, const double x) 
{
  int i, n_max;
  double r, v1;

  n_max = static_cast<int>(ceil(v - x - 1));
  v1 = v - 1.0;

  r = 1.0 / v1;
  for (i = 1; i <= n_max; i++) {
    r *= x / (v1 - i);
    if (r < epsilon)
      return i;
  }
  return 0;
}

double expint_asymp_v(const double v, const double x, const int n) 
{
  int i;
  double d, s, u, v1;

  v1 = v - 1.0;

  // series
  u = 1.0 / v1;
  s = u;
  d = 1.0;

  for (i = 1; i <= n; i++) {
    u /= (v1 - i);
    d *= -x;
    s += u * d;
  }
  return s * std::exp(-x);
}

// choose n such that the bound of the remainder of the series is less than
// the threshold = 2^(-53). Perform linear search (to be optimised).
int expint_series1_terms(const double v, const double x, double* bound, int* n) 
{
  int n_max = 50;
  double b, t, u, w;

  u = 1.0;
  t = 1.0 - v;
  w = 1.0;

  for (int i = 1; i < n_max; i++) {
    u *= x;
    t += 1;
    w *= i;
    b = u / (t * w);

    if (std::abs(b) < epsilon) {
      *n = i;
      *bound = b;
      return 1;
    }
  }
  return 0;
}

// choose n for series 2. Perform linear search
int expint_series2_terms(const double v, const double x, double* bound, int* n) 
{
  int n_max = 50;
  double b, t, u, w;

  // first iteration
  u = 1.0;
  t = 1.0 - v;
  w = 1.0;

  for (int i = 1; i < n_max; i++) {
    u *= x;
    w *= (t + i);
    b = u / w;

    if (std::abs(b) < epsilon) {
      *n = i;
      *bound = b;
      return 1;
    }
  }
  return 0;
}

// series expansion 1 for small x, alternating series.
double expint_series_a(const double v, const double x) 
{
  int k, n = 0;
  double b, con, d, q, r, u, shi, slo, t, thi, tlo;

  // select number of terms
  expint_series1_terms(v, x, &b, &n);

  t = 1.0 - v;
  u = std::tgamma(t) * std::pow(x, -t);

  shi = 1.0 / t;
  slo = 0.0;
  q = 1.0;
  d = 1.0;

  for (k = 1; k <= n; k++) {
    q *= -x;
    d *= k;
    r = q / (d * (t + k));

    // sum double-double
    thi = shi + r;
    con = thi - shi;
    tlo = (shi - (thi - con) + (r - con));
    tlo += slo;

    shi = thi + tlo;
    slo = (thi - shi) + tlo;
  }
  shi = -shi;
  slo = -slo;
  thi = shi + u;
  con = thi - shi;
  tlo = (shi - (thi - con) + (u - con));
  tlo += slo;

  return thi + tlo;
}

double expint_series_b(const double v, const double x) 
{
  int k, n = 0;
  double aux, b, con, d, q, r, u, shi, slo, t, thi, tlo;

  // select number of terms
  expint_series2_terms(v, x, &b, &n);

  t = 1.0 - v;
  u = std::tgamma(t) * std::pow(x, -t);

  shi = 1.0;
  slo = 0.0;
  q = 1.0;
  d = 1.0;

  for (k = 1; k <= n; k++) {
    q *= x;
    d *= (t + k);
    r = q / d;
    // sum double-double
    thi = shi + r;
    con = thi - shi;
    tlo = (shi - (thi - con) + (r - con));
    tlo += slo;

    shi = thi + tlo;
    slo = (thi - shi) + tlo;
  }
  aux = exp(-x) / (v - 1.0);
  shi = shi * aux;
  slo = slo * aux;
  thi = shi + u;
  con = thi - shi;
  tlo = (shi - (thi - con) + (u - con));
  tlo += slo;

  return thi + tlo;
}

// compute E_n(x), n is integer positive
// code based on cephes expn.c power series expansion
double expint_series_n(const int n, const double x)
{
  int i, k, terms = 0;
  double psi0, factor;
  double b, xk, yk, pk;
  double r, shi, slo, thi, tlo, con;

  // select number of terms
  expint_series1_terms(static_cast<double>(n), x, &b, &terms);

  // compute digamma function, \psi_0, using its series expansion for integer
  // argument:
  // \psi_0 = -EULER_MAS - sum_{i=1}^{n-1} 1/i
  psi0 = -GNSTLIB::constants::eulmasc;
  for (i = 1; i < n; i++)
    psi0 += 1.0 / i;

  // this series is used for n < 20, therefore a direct evaluation of the
  // factor is safe: (-x) ^ (n-1) * psi / gamma(n)
  factor = std::pow(-x, n - 1.0) * (psi0 - std::log(x)) / std::tgamma(n);

  // series
  xk = 0.0;
  yk = 1.0;
  pk = 1.0 - n;

  if (n == 1) {
    shi = 0.0;
    slo = 0.0;
  } 
  else {
    shi = 1.0 / pk;
    slo = 0.0;
  }

  for (k = 0; k <= terms; k++) {
    xk += 1.0;
    yk *= -x / xk;
    pk += 1.0;
    if (pk != 0.0) {
      r = yk / pk;
      thi = shi + r;
      // sum double-double
      con = thi - shi;
      tlo = (shi - (thi - con) + (r - con));
      tlo += slo;

      shi = thi + tlo;
      slo = (thi - shi) + tlo;
    }
  }
  shi = -shi;
  slo = -slo;
  thi = shi + factor;
  con = thi - shi;
  tlo = (shi - (thi - con) + (factor - con));
  tlo += slo;

  return thi + tlo;
}

// Convergent Laguerre series
double expint_laguerre_series(const double v, const double x) 
{
  const int maxiter = 500;
  const double tol = 0.1 * epsilon;

  int k;
  double Lk, Lk1, u, d, q, r, s;
  
  Lk = 1.0;
  Lk1 = x + v;
  // iteration 0
  s = 1/ Lk1;

  u = 1.0;
  d = 1.0;
  for (k = 1; k < maxiter; k++) {
    u *= v + k - 1;
    d *= 1 + k;
    q = (x + 2*k + v) / (k + 1) * Lk1 - (k + v - 1) / (k + 1) * Lk;

    r = u / (d * (q * Lk1));
    s += r;

    Lk = Lk1;
    Lk1 = q;
    if (std::abs(r) < tol)
        break;
  }
  return s * std::exp(-x);
}

// recursion starting with E1(x)
double expint_series_n_x_le_2(const int n, const double x) 
{
  int err_id, k, n1, m;
  double d, s, u;

  // compute E1(x)
  n1 = n - 1;

  // this series is used for n < 10, therefore a direct evaluation of
  // u is safe: (-x) ^ (n-1) * E1(x) / gamma(n)
  u = std::pow(-x, n1) * GNSTLIB::e1(x, err_id) / std::tgamma(n);

  // series
  m = n1;
  s = 1.0 / m;
  d = 1.0;
  for (k = 1; k <= n - 2; k++) {
      m *= n1 - k;
      d *= -x;
      s += d / m;
  }
  return u + std::exp(-x) * s;
}

double expint_n_e_acc(const int n, const double e, const double x)
{
  int k, terms = 0;
  double b;
   __float128 t, u, s, ss, q, d, r;
   __float128 one = 1.0q;

  expint_series1_terms(n+e, x, &b, &terms);

  t = 1 - n;
  t -= e;
  u = expq(lgammaq(t)); //tgammaq returns incorrect sign for very small e
  if (n % 2 == 0) {
    if (e < 0.0q)
      u = -fabsq(u);  
  }
  else {
    if (e > 0.0q)
      u = -fabsq(u);
  }

  u *= powq(x, -t);

  s = one / t;
  q = d = one;
  for (k = 1; k <= terms; k++)
  {
    q *= -x;
    d *= k;
    r = q / (d * (t + k));
    s += r;
  }
  ss = u - s;
  return static_cast<double>(ss);
}

// series expansion for small x
double expint_small_x(const double v, const double x) 
{
  if (v / x > 10.0) // fast convergence
    return expint_series_a(v, x);

  if (v > 1.5 && x > 0.5) // slow but accurate
    return expint_laguerre_series(v, x);
  else {
    if (v < 0.9) // all terms of the expansion are positive
      return expint_series_b(v, x);
    else
      return expint_series_a(v, x);
  }
}

// asymptotic expansion or Laguerre series for large x
double expint_large_x(const double v, const double x) 
{
  int iter;

  // fixed v and large v
  if (x / v > 100.0) {
    iter = use_expint_asymp_x(v, x);
    if (iter)
      return expint_asymp_x(v, x, iter);
    else
      return expint_laguerre_series(v, x);
  } 
  else
    return expint_laguerre_series(v, x);
}

// asympotitc expansion or Laguerre series for large v
double expint_large_v(const double v, const double x) 
{
  int iter;

  // fixed and small x and large v
  if (x < 5.0) {
    iter = use_expint_asymp_v(v, x);
    if (iter)
      return expint_asymp_v(v, x, iter);
  } 
  return expint_laguerre_series(v, x);
}

double expint_n(const int n, const double x)
{
  int err_id;

  // special cases
  if (n == 0)
    return std::exp(-x) / x;

  if (n == 1) {
    if (x > 0.9 && x < 10.0)
      return GNSTLIB::e1(x, err_id);
  }

  // small x
  if (x <= 1.5 && n < 20)
    return expint_series_n(n, x);
  else if (x <= 2.0 && n <= 10)
    // use series expansion in terms of E1(x) to avoid many iterations
    // laguerre series
    return expint_series_n_x_le_2(n, x);
  else if (n >= x)
    return expint_large_v(n, x);
  else
    return expint_large_x(n, x);
}

// Generalized exponential integral
double GNSTLIB::expint(const double v, const double x, int& err_id)
{
  const int vint = int(v);
  const double sqrtx = std::sqrt(x);

  double ans;

  err_id = 0;

  // error_id's
  if (v < 0) {
    err_id = 4;
    return nan;
  }

  if (x < 0) {
    err_id = 5;
    return nan;
  }

  // special cases
  if (x == 0.0) {
    if (v <= 1.0)
      return inf;
    else
      return 1.0 / (v - 1.0);
  }

  if (v == 0.5 && x < 10.0)
    return GNSTLIB::constants::sqrtpi * std::erfc(sqrtx) / sqrtx;

  // integer case
  if (v == vint) 
    ans = expint_n(vint, x);
  else {
    if (x <= 1.0)
      ans = expint_small_x(v, x);
    else if (v >= x)
      ans = expint_large_v(v, x);
    else
      ans = expint_large_x(v, x);
  }

  // floating-point error handler
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

double GNSTLIB::expint(const int n, const double e, const double x, int& err_id)
{
  double ans;

  // err_id's 
  if (n + e < 0) {
    err_id = 4;
    return nan;
  }

  if (std::abs(e) > 0.5) {
    err_id = 6;
    return nan;
  }

  // check if input is integer
  if (e == 0.0)
    ans = expint_n(n, x);

  // special case - higher precision is required
  if (std::abs(e) < 0.1 && n < 10 && x < 2.0)
    ans = expint_n_e_acc(n, e, x);
  else
    return GNSTLIB::expint(n + e, x, err_id);

  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}