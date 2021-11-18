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

#include <iomanip>
#include <iostream>

#include <gnstlib.hpp>

void basic_functions()
{
  int err = 0;
  
  std::cout << "\n-----------------------Basic Functions-----------------------" 
    << std::scientific << std::endl;
  std::cout << "expm1x(0.01)        = " << GNSTLIB::expm1x(0.01) << std::endl;
  std::cout << "log1pexp(3.2)       = " << GNSTLIB::log1pexp(3.2) << std::endl;
  std::cout << "log1mexp(-1.5)      = " << GNSTLIB::log1mexp(-1.5, err)
    << "\terr_id = " << err << std::endl;
  std::cout << "powm1(0.9, 0.3)     = " << GNSTLIB::powm1(0.9, 0.3, err) 
    << "\terr_id = " << err << std::endl;
}

void gamma_functions()
{
  int err = 0;
  double x, y;
  std::complex<double> z;

  x = 213.2; 
  y = 216.9;

  std::cout << "\n-----------------------Gamma Functions-----------------------"
    << std::scientific << std::endl;

  std::cout.precision(16);

  std::cout << "gamma(16.4)         = " << GNSTLIB::gamma(16.4, err);
  std::cout << "\terr_id = " << err << std::endl;
  std::cout << "qgamma(213.2,216.9) = " << GNSTLIB::qgamma(x, y, err)
    << "\terr_id = " << err << std::endl;
  std::cout << "gammastar(1431.2)   = " << GNSTLIB::gammastar(1431.2, err)
    << "\terr_id = " << err << std::endl;
  std::cout << "stirling(131.9)     = " << GNSTLIB::stirling(131.9, err)
    << "\terr_id = " << err << std::endl;

  z = std::complex<double>(0.0, 16.4);
  std::cout << "gamma(16.4i)        = " << std::setprecision(4) << 
    GNSTLIB::gamma(z, err) << "\terr_id = " << err << std::endl;
  z = std::complex<double>(6.0, 15.0);
  std::cout << "loggamma(6+15i)     = " << std::setprecision(4) << 
    GNSTLIB::loggamma(z, err) << "\terr_id = " << err << std::endl;
}

void exponential_integral_functions()
{
  int err = 0;
  std::complex<double> ci, si, z;

  z = std::complex<double>(2.0, 10.0);

  std::cout << "\n--------------------Exponential Integrals--------------------"
    << std::scientific << std::endl;

  std::cout.precision(16);

  std::cout << "ei(2.2)             = " << GNSTLIB::ei(2.2, err) 
    << "\terr_id = " << err << std::endl;
  std::cout << "e1(2.2)             = " << GNSTLIB::e1(2.2, err) 
    << "\terr_id = " << err << std::endl;
  std::cout << "li(200000)          = " << GNSTLIB::li(200000, err) 
    << "\terr_id = " << err << std::endl;

  // compute sine and cosine integral
  
  GNSTLIB::sici(z, si, ci, err);

  std::cout << "ci(2+10i)           = " << std::setprecision(4) << ci 
    << "\terr_id = " << err << std::endl;
  std::cout << "si(2+10i)           = " << std::setprecision(4) << si 
    << "\terr_id = " << err << std::endl;
}

void error_functions()
{
  int err = 0;
  std::complex<double> z;

  z = std::complex<double>(3.2, 6.5);

  std::cout << "\n-----------------------Error Functions-----------------------"
    << std::scientific << std::endl;

  std::cout.precision(16);

  std::cout << "erf(3.7)            = " << GNSTLIB::erf(3.7, err)
    << "\terr_id = " << err << std::endl;
  std::cout << "erfc(3.7)           = " << GNSTLIB::erfc(3.7, err)
    << "\terr_id = " << err << std::endl;
  std::cout << "fresnelc(55,9)      = " << GNSTLIB::fresnelc(55.9, err)
    << "\terr_id = " << err << std::endl;
  std::cout << "fresnels(55,9)      = " << GNSTLIB::fresnels(55.9, err)
    << "\terr_id = " << err << std::endl;

  std::cout << "faddeeva(3.2+6.5i)  = " << std::setprecision(4) << 
    GNSTLIB::faddeeva(z, err) << "\terr_id = " << err << std::endl;  
}

void incomplete_gamma_functions()
{
  int err = 0;

  std::cout << "\n----Incomplete Gamma and Generalized Exponential Integral----"
    << std::scientific << std::endl;

  std::cout.precision(16);

  std::cout << "gammainc_p(2.3,8.3) = " << GNSTLIB::gammainc_p(2.3, 8.3, err)
    << "\terr_id = " << err << std::endl;

  std::cout << "gammainc_q(2.3,8.3) = " << GNSTLIB::gammainc_q(2.3, 8.3, err)
    << "\terr_id = " << err << std::endl;

  std::cout << "expint(120.2, 12.2) = " << GNSTLIB::expint(120.2, 12.2, err)
    << "\terr_id = " << err << std::endl;
}

void vectorized_functions()
{
  std::cout << "\n-------------Example vectorized functions: Gamma-------------"
    << std::scientific << std::endl;

  std::cout.precision(16);

  const int n = 3;

  std::vector<double> v = {1.0, 2.3, 4.3};
  std::vector<double> r;

  GNSTLIB::gamma_vec(v, r, 0);

  for (auto& i:r) {
    std::cout << "gamma(" << v[&i-&r[0]] << ") = " << i << std::endl;
  }

  const double varray[] = {1.0, 2.3, 4.3};
  double rarray[n];

  GNSTLIB::gamma_vec(n, varray, rarray, 1);

  for (int i = 0; i < n; i++)
    std::cout << "gamma(" << varray[i] << ") = " << rarray[i] << std::endl;
}

int main()
{
  std::cout.precision(16);

  std::cout << "GNSTLIB version 0.1 - Release October 2017" << std::endl;
  std::cout << "==========================================" << std::endl;


  basic_functions();
  gamma_functions();
  exponential_integral_functions();
  error_functions();
  incomplete_gamma_functions();
  vectorized_functions();
}