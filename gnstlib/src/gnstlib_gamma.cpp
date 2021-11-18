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

#include "gnstlib.hpp"
#include "gnstlib_polyeval.hpp"
#include "gnstlib_utils.hpp"

typedef std::complex<double> cmplx;

// Gamma function for real argument
double GNSTLIB::gamma(const double x, int& err_id) 
{
  double ans;

  err_id = 0;

  // special cases
  if (x > 0.0 && x < 1.0 / maxdouble) {
    // overflow
    err_id = 1;
    return inf;
  } else if (x <= 0.0 && x == std::floor(x)) {
    // zero, negative integer or -inf -> invalid argument
    err_id = 4;
    return inf;
  }

  ans = std::tgamma(x);

  // floating-point error handler
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}


// Logarithm of the absolute value of the Gamma function for real argument
double GNSTLIB::gammaln(const double x, int& err_id) 
{
  double ans;

  err_id = 0;  

  // special cases
  if (x <= 0.0 && x == std::floor(x)) {
    // zero, negative integer or -inf -> invalid argument
    err_id = 4;
    return inf;
  }

  ans = std::lgamma(x);

  // floating-point error handler
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}


// Sign of the Gamma function - based on Scipy implementation
int GNSTLIB::gammasign(const double x) 
{
  double fx;

  fx = std::floor(x);

  if (x > 0.0)
    return 1;
  else
    if (x - fx == 0.0)
      return 0;
    else if (static_cast<int>(fx) % 2)
      return -1;
    else
      return 1;
}

cmplx loggamma_asymptotic(const cmplx z) 
{
  cmplx poly, rz, rzz;
  const double hlog2pi = 0.918938533204672742;

  std::vector<double> coeffs = {
    8.3333333333333333333e-2,  //  1/12
    -2.7777777777777777778e-3, // -1/360
    7.9365079365079365079e-4,  //  1/1260
    -5.952380952380952381e-4,  // -1/1680
    8.4175084175084175084e-4,  //  1/1188
    -1.9175269175269175269e-3, // -691/360360
    6.4102564102564102564e-3,  //  1/156
    -2.955065359477124183e-2,  // -3617/122400
  };

  rz = 1.0 / z;
  rzz = rz / z;

  return (z - 0.5) * std::log(z) - z + hlog2pi + rz * polyeval(coeffs, rzz);
}

cmplx loggamma_taylor(const cmplx z) 
{
  std::vector<double> coeffs = {
    -5.7721566490153286061e-1, 8.2246703342411321824e-1,
    -4.0068563438653142847e-1, 2.7058080842778454788e-1,
    -2.0738555102867398527e-1, 1.6955717699740818995e-1,
    -1.4404989676884611812e-1, 1.2550966952474304242e-1,
    -1.1133426586956469049e-1, 1.0009945751278180853e-1,
    -9.0954017145829042233e-2, 8.3353840546109004025e-2,
    -7.6932516411352191473e-2, 7.1432946295361336059e-2,
    -6.6668705882420468033e-2, 6.2500955141213040742e-2,
    -5.8823978658684582339e-2, 5.5555767627403611102e-2,
    -5.2631679379616660734e-2, 5.000004769810169364e-2,
    -4.7619070330142227991e-2, 4.5454556293204669442e-2,
    -4.3478266053040259361e-2
  };

  cmplx u = z - 1.0;
  return u * polyeval(coeffs, u);
}

cmplx loggamma_recurrence(const cmplx z) 
{
  int signflips = 0;
  int sb = 0;
  int nsb;
  cmplx ans;
  cmplx shiftprod = z;
  double x = real(z);

  x += 1.0;
  while (x <= 7.0) {
    shiftprod *= cmplx(x, imag(z));
    nsb = std::signbit(imag(shiftprod));
    signflips += (nsb != 0 && sb == 0) ? 1 : 0;
    sb = nsb;
    x += 1.0;
  }
  ans = loggamma_asymptotic(cmplx(x, imag(z))) - std::log(shiftprod);
  ans -= cmplx(0.0, signflips * GNSTLIB::constants::twopi);
  return ans;
}

/**
 * Principal branch of the logarithm of the Gamma function
 * 
 * References:
 *  https://github.com/JuliaLang/julia/blob/master/base/special/gamma.jl
 *  https://github.com/scipy/scipy/special/_loggamma.pxd 
 */
cmplx GNSTLIB::loggamma(const cmplx z, int& err_id) 
{
  double re_z = real(z);
  double im_z = imag(z);
  double zaux0, zaux1;

  err_id = 0;

  // special cases
  if (std::isnan(re_z) || std::isnan(im_z)) {
    err_id = -1;
    return cmplx(nan, nan);
  }
  else if (re_z <= 0.0 && z == std::floor(re_z)) {
    err_id = 4;
    return cmplx(nan, nan);
  }
  else if (re_z > 0.0 && imag(z) == 0.0)
    return GNSTLIB::gammaln(re_z, err_id);
  else if (re_z > 7.0 || std::abs(im_z) > 7.0){
    return loggamma_asymptotic(z);
  }
  else if (std::abs(z - 1.0) <= 0.2){
    return loggamma_taylor(z);
  }
  else if (std::abs(z - 2.0) <= 0.2) {
    return std::log(z - 1.0) + loggamma_taylor(z - 1.0);
  }
  else if (re_z < 0.1) {
    // reflection formula
    zaux0 = std::floor(0.5 * re_z + 0.25);
    zaux1 = std::copysign(GNSTLIB::constants::twopi, im_z) * zaux0;
    return cmplx(GNSTLIB::constants::lnpi, zaux1) 
      - std::log(GNSTLIB::sinpi(z)) - GNSTLIB::loggamma(1.0 - z, err_id);
  } else if (std::signbit(im_z) == 0){
    return loggamma_recurrence(z);
  }
  else{
    return std::conj(loggamma_recurrence(cmplx(re_z, -im_z)));
  }
}

// Gamma function for complex argument
cmplx GNSTLIB::gamma(const cmplx z, int& err_id) 
{
  err_id = 0;

  if (real(z) <= 0.0 && z == std::floor(real(z))) {
    err_id = 4;
    return cmplx(nan, nan);
  }
  return std::exp(GNSTLIB::loggamma(z, err_id));
}


 // Reciprocal of Gamma function
cmplx GNSTLIB::rgamma(const cmplx z) 
{
  int err_id = 0; // dummy argument

  if (real(z) <= 0.0 && z == std::floor(real(z)))
    return 0.0;
  
  return std::exp(-GNSTLIB::loggamma(z, err_id));
}

// axugam(x) = 1 / gamma(1 + x), x \in [-1, 1]
double GNSTLIB::auxgam(const double x, int& err_id)
{
  double ans, t, xp1;

  if (std::abs(x) > 1.0) {
    err_id = 4;
    return nan;
  }

  xp1 = 1.0 + x;

  if (x < 0.0)
    ans = -(1.0 + xp1 * xp1 * GNSTLIB::auxgam(xp1, err_id)) / (1.0 - x);
  else {
    std::vector<double> dr = {
      -1.013609258009865776949,
      0.784903531024782283535e-1,
      0.67588668743258315530e-2,
      -0.12790434869623468120e-2,
      0.462939838642739585e-4,
      0.43381681744740352e-5,
      -0.5326872422618006e-6,
      0.172233457410539e-7,
      0.8300542107118e-9,
      -0.10553994239968e-9,
      0.39415842851e-11,
      0.362068537e-13,
      -0.107440229e-13,
      0.5000413e-15,
      -0.62452e-17,
      -0.5185e-18,
      0.347e-19,
      -0.9e-21      
    };
    t = 2.0 * x - 1.0;
    ans = chepolsum(dr, t);
  }
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// log(gamma(1 + x)) for x \in [-1, 1]
double lngam1(const double x)
{
  int err_id;

  return -std::log1p(x * (x - 1.0) * GNSTLIB::auxgam(x, err_id));
}

// Stirling series for x > 0
double GNSTLIB::stirling(const double x, int& err_id)
{
  double ans, c0, c1, c2, c3, c4, c5, c6, z;

  err_id = 0;

  // special cases
  if (x > 0.0 && x < 1.0 / maxdouble) {
    err_id = 1;
    return inf;
  } else if (x <= 0.0) {
    err_id = 4;
    return nan;
  }

  else if (x < 1.0) 
    ans = lngam1(x) - (x + 0.5) * std::log(x) + x 
        - GNSTLIB::constants::lnsqrttwopi;
  else if (x < 2.0)
    ans = lngam1(x - 1.0) - (x - 0.5) * std::log(x) + x 
        - GNSTLIB::constants::lnsqrttwopi;
  else if (x < 3.0)
    ans = lngam1(x - 2.0) - (x - 0.5) * std::log(x) + x 
        - GNSTLIB::constants::lnsqrttwopi + std::log(x - 1.0);
  else if (x < 12.0) {
    std::vector<double> a = {
      1.996379051590076518221,
      -0.17971032528832887213e-2,
      0.131292857963846713e-4,
      -0.2340875228178749e-6,
      0.72291210671127e-8,
      -0.3280997607821e-9,
      0.198750709010e-10,
      -0.15092141830e-11,
      0.1375340084e-12,
      -0.145728923e-13,
      0.17532367e-14,
      -0.2351465e-15,
      0.346551e-16,
      -0.55471e-17,
      0.9548e-18,
      -0.1748e-18,
      0.332e-19,
      -0.58e-20
    };

    z = 18.0 / (x * x) - 1.0;
    ans = chepolsum(a, z) / (12.0 * x);
  } else {
    z = 1.0 / (x * x);
    if (x < 1000.0) {
      c0=0.25721014990011306473e-1;
      c1=0.82475966166999631057e-1;
      c2=-0.25328157302663562668e-2;
      c3=0.60992926669463371e-3;
      c4=-0.33543297638406e-3;
      c5=0.250505279903e-3;
      c6=0.30865217988013567769; 

      ans = ((((((c5*z+c4)*z+c3)*z+c2)*z+c1)*z+c0)/(c6+z)/x);

    } else {
      ans = (((-z * 0.000595238095238095238095238095238 
        + 0.000793650793650793650793650793651) * z
        - 0.00277777777777777777777777777778) * z
        + 0.0833333333333333333333333333333) / x;
    }
  }

  // floating-point error handler
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// Computation of the scaled gamma function for x > 0
double GNSTLIB::gammastar(const double x, int& err_id)
{
  double ans, t;

  if (x <= 0) {
    err_id = 4;
    return nan;
  }
  if (x >= 3.0)
    ans = std::exp(GNSTLIB::stirling(x, err_id));
  else {
    t = std::exp(-x + (x + 0.5) * std::log(x)) * GNSTLIB::constants::sqrttwopi;
    ans = GNSTLIB::gamma(x, err_id) / t;
  }

  // floating-point error handler
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;
}

// Computation of the quotient of two gamma functions gamma(x) / gamma(y)
// requires: qratio, shiftfact, qasym

double qratio(const double x, const double y)
{
  int err_id;
  double b, c, q, t, z1, z2;

  if ((x <= 0.0 || y <= 0.0) || y > 1.5 * x)
    if (x < -170.0 || y < -170.0) {
      z1 = 1.0 - x;
      z2 = 1.0 - y;
      t = GNSTLIB::sinpi(z2) / GNSTLIB::sinpi(z1);
      return t * GNSTLIB::qgamma(z2, z1, err_id);
    }
    else
      return GNSTLIB::gamma(x, err_id) / GNSTLIB::gamma(y, err_id);
  else {
    b = y - x;
    c = b / x;
    q = GNSTLIB::log1pmx(c);
    q = std::exp(-b * std::log(x) - (b - 0.5) * std::log(1.0 + c) - x * q);
    return q * GNSTLIB::gammastar(x, err_id) / GNSTLIB::gammastar(y, err_id);
  }
}

double qasym(const double x, const double y)
{
  int k;
  double r, r2, r3, r4, r5, r6, r7, s, t, u, w, w2;

  w = 0.5 * (x + y - 1.0);
  w2 = 1.0 / (w * w);
  r = 0.5 * (x - y + 1.0);
  r2 = r * r; r3 = r * r2; r4 = r * r3; r5 = r * r4; r6 = r * r5; r7 = r * r6;

  std::vector<double> cc = {
    1.0,
    8.333333333333333333e-02 * r,
    6.944444444444444444e-04 * r  + 3.472222222222222222e-03 * r2,
    1.102292768959435626e-05 * r  + 5.787037037037037037e-05 * r2 + 
    9.645061728395061728e-05 * r3,
    2.066798941798941799e-07 * r  + 1.159703850676072898e-06 * r2 + 
    2.411265432098765432e-06 * r3 + 2.009387860082304527e-06 * r4,
    4.175351397573619796e-09 * r  + 2.487813541054281795e-08 * r2 + 
    5.836793307858122673e-08 * r3 + 6.697959533607681756e-08 * r4 + 
    3.348979766803840878e-08 * r5,
    8.806983564479155308e-11 * r  + 5.522261215144078166e-10 * r2 + 
    1.411355758867332941e-09 * r3 + 1.900413121638687482e-09 * r4 + 
    1.395408236168267032e-09 * r5 + 4.651360787227556775e-10 * r6,
    1.911790932954954119e-12 * r  + 1.251692008048563016e-11 * r2 + 
    3.417901600113656747e-11 * r3 + 5.116496865950312452e-11 * r4 +
    4.540614101817376852e-11 * r5 + 2.325680393613778387e-11 * r6 + 
    5.537334270508996161e-12 * r7,
  };

  s = 1.0;
  t = 1.0;
  u = 1.0;

  k = 1;
  do {
    t = -4.0 * w2 * (k - r) * (k - r - 0.5) * t;
    u = t * cc[k];
    s += u;
    k += 1;
  } while (std::abs(u) > epsilon && k < 8);

  return s * std::exp((x - y) * std::log(w));
}

double shiftfact(const double x, const int n)
{
  int k;
  double s;

  if (n == 0)
    return 1.0;
  else if (n < 0)
    return 1.0 / shiftfact(x - n, n);
  else {
    s = 1.0;
    for (k = 0; k < n; k++) 
      s *= (x + k);
    return s;
  }
}

double GNSTLIB::qgamma(const double x, const double y, int& err_id)
{
  double ans;
  int n;

  if (x <= 0.0 || y <= 0.0) {
    if (x == std::floor(x) || y == std::floor(y)) {
      err_id = 4;
      return nan;
    }
    ans = qratio(x, y);
  } else if (x > y)
    ans = 1.0 / GNSTLIB::qgamma(y, x, err_id);
  else {
    n = static_cast<int>(y - x);
    if (n > 15)
      if ((std::floor(x) == x && std::floor(y) == y) && (x <= 50 && y <= 50)) {
        ans = GNSTLIB::factorial(static_cast<int>(x)-1, err_id) / 
              GNSTLIB::factorial(static_cast<int>(y)-1, err_id);
      }
      else
        ans = qratio(x, y);
    else if (n > 0) {
      if (n == (y - x) && n < 10)
        ans = 1.0 / shiftfact(x, n); // faster, loss of accuracy as n grows
      else
        ans = GNSTLIB::qgamma(x + n, y, err_id) / shiftfact(x, n);
    }
    else if (x < 26.0)
      ans = qratio(x, y);
    else{
      ans = qasym(x, y);
    }
  }

  // floating-point error handler
  err_id = gnstlib_fp_error_handler(ans); 
  return ans;  
}

// factorial(n)
double GNSTLIB::factorial(const int n, int& err_id) 
{  
  err_id = 0;

  // special cases
  if (n > 170){
    err_id = 1; 
    return inf;
  }

  if (n < 0.0) {
    err_id = 4;
    return nan;
  }

  const std::vector<double> factorials = {
      1,
      1,
      2,
      6,
      24,
      120,
      720,
      5040,
      40320,
      362880.0L,
      3628800.0L,
      39916800.0L,
      479001600.0L,
      6227020800.0L,
      87178291200.0L,
      1307674368000.0L,
      20922789888000.0L,
      355687428096000.0L,
      6402373705728000.0L,
      121645100408832000.0L,
      0.243290200817664e19L,
      0.5109094217170944e20L,
      0.112400072777760768e22L,
      0.2585201673888497664e23L,
      0.62044840173323943936e24L,
      0.15511210043330985984e26L,
      0.403291461126605635584e27L,
      0.10888869450418352160768e29L,
      0.304888344611713860501504e30L,
      0.8841761993739701954543616e31L,
      0.26525285981219105863630848e33L,
      0.822283865417792281772556288e34L,
      0.26313083693369353016721801216e36L,
      0.868331761881188649551819440128e37L,
      0.29523279903960414084761860964352e39L,
      0.103331479663861449296666513375232e41L,
      0.3719933267899012174679994481508352e42L,
      0.137637530912263450463159795815809024e44L,
      0.5230226174666011117600072241000742912e45L,
      0.203978820811974433586402817399028973568e47L,
      0.815915283247897734345611269596115894272e48L,
      0.3345252661316380710817006205344075166515e50L,
      0.1405006117752879898543142606244511569936e52L,
      0.6041526306337383563735513206851399750726e53L,
      0.265827157478844876804362581101461589032e55L,
      0.1196222208654801945619631614956577150644e57L,
      0.5502622159812088949850305428800254892962e58L,
      0.2586232415111681806429643551536119799692e60L,
      0.1241391559253607267086228904737337503852e62L,
      0.6082818640342675608722521633212953768876e63L,
      0.3041409320171337804361260816606476884438e65L,
      0.1551118753287382280224243016469303211063e67L,
      0.8065817517094387857166063685640376697529e68L,
      0.427488328406002556429801375338939964969e70L,
      0.2308436973392413804720927426830275810833e72L,
      0.1269640335365827592596510084756651695958e74L,
      0.7109985878048634518540456474637249497365e75L,
      0.4052691950487721675568060190543232213498e77L,
      0.2350561331282878571829474910515074683829e79L,
      0.1386831185456898357379390197203894063459e81L,
      0.8320987112741390144276341183223364380754e82L,
      0.507580213877224798800856812176625227226e84L,
      0.3146997326038793752565312235495076408801e86L,
      0.1982608315404440064116146708361898137545e88L,
      0.1268869321858841641034333893351614808029e90L,
      0.8247650592082470666723170306785496252186e91L,
      0.5443449390774430640037292402478427526443e93L,
      0.3647111091818868528824985909660546442717e95L,
      0.2480035542436830599600990418569171581047e97L,
      0.1711224524281413113724683388812728390923e99L,
      0.1197857166996989179607278372168909873646e101L,
      0.8504785885678623175211676442399260102886e102L,
      0.6123445837688608686152407038527467274078e104L,
      0.4470115461512684340891257138125051110077e106L,
      0.3307885441519386412259530282212537821457e108L,
      0.2480914081139539809194647711659403366093e110L,
      0.188549470166605025498793226086114655823e112L,
      0.1451830920282858696340707840863082849837e114L,
      0.1132428117820629783145752115873204622873e116L,
      0.8946182130782975286851441715398316520698e117L,
      0.7156945704626380229481153372318653216558e119L,
      0.5797126020747367985879734231578109105412e121L,
      0.4753643337012841748421382069894049466438e123L,
      0.3945523969720658651189747118012061057144e125L,
      0.3314240134565353266999387579130131288001e127L,
      0.2817104114380550276949479442260611594801e129L,
      0.2422709538367273238176552320344125971528e131L,
      0.210775729837952771721360051869938959523e133L,
      0.1854826422573984391147968456455462843802e135L,
      0.1650795516090846108121691926245361930984e137L,
      0.1485715964481761497309522733620825737886e139L,
      0.1352001527678402962551665687594951421476e141L,
      0.1243841405464130725547532432587355307758e143L,
      0.1156772507081641574759205162306240436215e145L,
      0.1087366156656743080273652852567866010042e147L,
      0.103299784882390592625997020993947270954e149L,
      0.9916779348709496892095714015418938011582e150L,
      0.9619275968248211985332842594956369871234e152L,
      0.942689044888324774562618574305724247381e154L,
      0.9332621544394415268169923885626670049072e156L,
      0.9332621544394415268169923885626670049072e158L,
      0.9425947759838359420851623124482936749562e160L,
      0.9614466715035126609268655586972595484554e162L,
      0.990290071648618040754671525458177334909e164L,
      0.1029901674514562762384858386476504428305e167L,
      0.1081396758240290900504101305800329649721e169L,
      0.1146280563734708354534347384148349428704e171L,
      0.1226520203196137939351751701038733888713e173L,
      0.132464181945182897449989183712183259981e175L,
      0.1443859583202493582204882102462797533793e177L,
      0.1588245541522742940425370312709077287172e179L,
      0.1762952551090244663872161047107075788761e181L,
      0.1974506857221074023536820372759924883413e183L,
      0.2231192748659813646596607021218715118256e185L,
      0.2543559733472187557120132004189335234812e187L,
      0.2925093693493015690688151804817735520034e189L,
      0.339310868445189820119825609358857320324e191L,
      0.396993716080872089540195962949863064779e193L,
      0.4684525849754290656574312362808384164393e195L,
      0.5574585761207605881323431711741977155627e197L,
      0.6689502913449127057588118054090372586753e199L,
      0.8094298525273443739681622845449350829971e201L,
      0.9875044200833601362411579871448208012564e203L,
      0.1214630436702532967576624324188129585545e206L,
      0.1506141741511140879795014161993280686076e208L,
      0.1882677176888926099743767702491600857595e210L,
      0.237217324288004688567714730513941708057e212L,
      0.3012660018457659544809977077527059692324e214L,
      0.3856204823625804217356770659234636406175e216L,
      0.4974504222477287440390234150412680963966e218L,
      0.6466855489220473672507304395536485253155e220L,
      0.8471580690878820510984568758152795681634e222L,
      0.1118248651196004307449963076076169029976e225L,
      0.1487270706090685728908450891181304809868e227L,
      0.1992942746161518876737324194182948445223e229L,
      0.269047270731805048359538766214698040105e231L,
      0.3659042881952548657689727220519893345429e233L,
      0.5012888748274991661034926292112253883237e235L,
      0.6917786472619488492228198283114910358867e237L,
      0.9615723196941089004197195613529725398826e239L,
      0.1346201247571752460587607385894161555836e242L,
      0.1898143759076170969428526414110767793728e244L,
      0.2695364137888162776588507508037290267094e246L,
      0.3854370717180072770521565736493325081944e248L,
      0.5550293832739304789551054660550388118e250L,
      0.80479260574719919448490292577980627711e252L,
      0.1174997204390910823947958271638517164581e255L,
      0.1727245890454638911203498659308620231933e257L,
      0.2556323917872865588581178015776757943262e259L,
      0.380892263763056972698595524350736933546e261L,
      0.571338395644585459047893286526105400319e263L,
      0.8627209774233240431623188626544191544816e265L,
      0.1311335885683452545606724671234717114812e268L,
      0.2006343905095682394778288746989117185662e270L,
      0.308976961384735088795856467036324046592e272L,
      0.4789142901463393876335775239063022722176e274L,
      0.7471062926282894447083809372938315446595e276L,
      0.1172956879426414428192158071551315525115e279L,
      0.1853271869493734796543609753051078529682e281L,
      0.2946702272495038326504339507351214862195e283L,
      0.4714723635992061322406943211761943779512e285L,
      0.7590705053947218729075178570936729485014e287L,
      0.1229694218739449434110178928491750176572e290L,
      0.2004401576545302577599591653441552787813e292L,
      0.3287218585534296227263330311644146572013e294L,
      0.5423910666131588774984495014212841843822e296L,
      0.9003691705778437366474261723593317460744e298L,
      0.1503616514864999040201201707840084015944e301L,
      0.2526075744973198387538018869171341146786e303L,
      0.4269068009004705274939251888899566538069e305L,
      0.7257415615307998967396728211129263114717e307L,
  };

  return factorials[n];
}
