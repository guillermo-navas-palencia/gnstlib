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
 * gnstlib_constants.hpp
 *  + machine constants declaration
 *  + mathematical constants
 */

#ifndef GNSTLIB_CONSTANTS_H
#define GNSTLIB_CONSTANTS_H

#include <limits>

// machine constants - all must be machine dependent

//machine precision
#define epsilon       std::numeric_limits<double>::epsilon()

// common limits
#define inf           std::numeric_limits<double>::infinity()
#define nan           std::numeric_limits<double>::quiet_NaN()
#define maxdouble     std::numeric_limits<double>::max()
#define mindouble     std::numeric_limits<double>::min()

// safe underflow and overflow limit
#define dwarf         std::numeric_limits<double>::min() * 1000.0
#define giant         std::numeric_limits<double>::max() / 1000.0 

// specific limits
#define explow       -300.0
#define exphigh       300.0


namespace GNSTLIB 
{

namespace constants 
{

// common constants
const double pi            = 3.1415926535897932385; // pi
const double pihalf        = 1.5707963267948966192; // pi/2
const double piquart       = 0.7853981633974483096; // pi/4
const double pithreequart  = 2.3561944901923449288; // 3pi/4
const double sqrtpi        = 1.7724538509055160272; // sqrt(pi)
const double lnpi          = 1.1447298858494001741; // log(pi)
const double twopi         = 6.2831853071795864769; // 2*pi
const double oneopi        = 0.3183098861837906715; // 1/pi
const double oneosqrtpi    = 0.5641895835477562869; // 1/sqrt(pi)
const double twoosqrtpi    = 1.1283791670955125739; // 2/sqrt(pi)
const double sqrttwopi     = 2.5066282746310005024; // sqrt(2*pi)
const double sqrttwoopi    = 0.7978845608028653559; // sqrt(2/pi)
const double lnsqrttwopi   = 0.9189385332046727418; // log(sqrt(2*pi))
const double onethird      = 0.3333333333333333333; // 1/3
const double twothird      = 0.6666666666666666667; // 2/3
const double onesix        = 0.1666666666666666667; // 1/6
const double twoexp14      = 1.1892071150027210667; // 2**(1/4)
const double twoexp13      = 1.2599210498948731648; // 2**(1/3)
const double sqrt2         = 1.4142135623730950488; // sqrt(2)
const double sqrt3         = 1.7320508075688772935; // sqrt(3)
const double log_2         = 0.6931471805599453094; // log(2)
const double log_10        = 2.3025850929940456840; // log(10)

// other constants and trascendental numbers
const double eulmasc       = 0.5772156649015328606; // Euler-Mascheroni
const double e             = 2.7182818284590452354; // e = exp(1) number
const double phi           = 1.6180339887498948482; // phi = 1/2(1+sqrt(5))
const double catalan       = 0.9159655941772190151; // Catalan constant
const double apery         = 1.2020569031595942854; // Apery constant
const double khinchin      = 2.6854520010653064453; // Khinchin constant
const double glaisher      = 1.2824271291006226369; // Glaisher constant
const double mertens       = 0.2614972128476427838; // Mertens constant
const double twinprime     = 0.6601618158468695739; // Twin primes constant
const double brun          = 1.902160583104;        // Brun's constant
const double mrb           = 0.1878596424620671202; // MRB constant 
const double madelung      = -1.747564594633182191; // Madelung constant
const double delian        = 1.2599210498948731648; // Delian constant
const double soldner       = 1.4513692348833810503; // Soldner constant

// Stieltjes numbers
const double stieltjesgamma0  = 0.5772156649015329;
const double stieltjesgamma1  = -0.07281584548367672;
const double stieltjesgamma2  = -0.009690363192872318;
const double stieltjesgamma3  = 0.002053834420303346;
const double stieltjesgamma4  = 0.002325370065467300;
const double stieltjesgamma5  = 0.0007933238173010627;
const double stieltjesgamma6  = -0.0002387693454301996;
const double stieltjesgamma7  = -0.0005272895670577510;
const double stieltjesgamma8  = -0.0003521233538030395;
const double stieltjesgamma9  = -0.00003439477441808805;
const double stieltjesgamma10 = 0.0002053328149090648;

} // namespace constants
} // namespace GNSTLIB

#endif /* GNSTLIB_CONSTANTS_H */
