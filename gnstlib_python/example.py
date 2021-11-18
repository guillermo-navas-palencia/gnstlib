# Copyright (C) 2017 A. Gil, G. Navas-Palencia, J. Segura and N. M. Temme
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
# 
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

import gnstlib


def basic_functions():
    print("\n-----------------------Basic Functions-----------------------")

    print("expm1x(0.01)         = {:.16e}".format(gnstlib.expm1x(0.01)))
    print("log1pexp(3.2)        = {:.16e}".format(gnstlib.log1pexp(3.2)))
    _log1mexp = gnstlib.log1mexp(-1.5)
    print("log1mexp(-1.5)       = {0:.16e}\t err_id = {1}".format(_log1mexp[0], 
        _log1mexp[1]))
    _powm1 = gnstlib.powm1(0.9, 0.3)
    print("powm1(0.9, 0.3)      = {0:.16e}\t err_id = {1}".format(_powm1[0], 
        _powm1[1]))


def gamma_functions():
    print("\n-----------------------Gamma Functions-----------------------")

    _gamma = gnstlib.gamma(16.4)
    print("gamma(16.4)          = {0:.16e}\t err_id = {1}".format(_gamma[0], 
        _gamma[1]))
    _qgamma = gnstlib.qgamma(213.2, 216.9)
    print("qgamma(213.2, 216.9) = {0:.16e}\t err_id = {1}".format(_qgamma[0], 
        _qgamma[1]))
    _gammastar = gnstlib.gammastar(1431.2)
    print("gammastar(1431.2)    = {0:.16e}\t err_id = {1}".format(_gammastar[0], 
        _gammastar[1]))
    _stirling = gnstlib.stirling(131.9)
    print("stirling(131.9)      = {0:.16e}\t err_id = {1}".format(_stirling[0], 
        _stirling[1]))

    _gamma = gnstlib.gamma(16.4j)
    re_gamma = _gamma[0].real
    im_gamma = _gamma[0].imag
    print("gamma(16.4j)         = {0:.16e}\n \
               {1:.16e}j\t err_id = {2}".format(re_gamma, im_gamma, _gamma[1]))

    _loggamma = gnstlib.loggamma(6+15j)
    re_loggamma = _loggamma[0].real
    im_loggamma = _loggamma[0].imag
    print("loggamma(6+15j)      = {0:.16e}\n \
               {1:.16e}j\t err_id = {2}".format(re_loggamma, im_loggamma, 
                _loggamma[1]))


def exponential_integral_functions():
    print("\n--------------------Exponential Integrals--------------------")

    _ei = gnstlib.ei(2.2)
    print("ei(2.2)              = {0:.16e}\t err_id = {1}".format(_ei[0], 
        _ei[1]))
    _e1 = gnstlib.e1(2.2)
    print("e1(2.2)              = {0:.16e}\t err_id = {1}".format(_e1[0], 
        _e1[1]))
    _li = gnstlib.li(200000)
    print("li(200000)           = {0:.16e}\t err_id = {1}".format(_li[0], 
        _li[1]))


def error_functions():
    print("\n-----------------------Error Functions-----------------------")

    _erf = gnstlib.erf(3.7)
    print("erf(2.2)             = {0:.16e}\t err_id = {1}".format(_erf[0], 
        _erf[1]))
    _erfc = gnstlib.erfc(3.7)
    print("erfc(2.2)            = {0:.16e}\t err_id = {1}".format(_erfc[0], 
        _erfc[1]))
    _fresnelc = gnstlib.fresnelc(55.9)
    print("fresnelc(55.9)       = {0:.16e}\t err_id = {1}".format(_fresnelc[0], 
        _fresnelc[1]))
    _fresnels = gnstlib.fresnels(55.9)
    print("fresnels(55.9)       = {0:.16e}\t err_id = {1}".format(_fresnels[0], 
        _fresnels[1]))
    _faddeeva = gnstlib.faddeeva(3.2+6.5j)
    re_faddeeva = _faddeeva[0].real
    im_faddeeva = _faddeeva[0].imag

    print("faddeeva(3.2+6.5j)   = {0:.16e}\n \
               {1:.16e}j\t err_id = {2}".format(re_faddeeva, im_faddeeva, 
                _faddeeva[1]))  


def incomplete_gamma_functions():
    print("\n----Incomplete Gamma and Generalized Exponential Integral----")

    _incp = gnstlib.gammainc_p(2.3, 8.3)
    print("gammainc_p(2.3, 8.3) = {0:.16e}\t err_id = {1}".format(_incp[0], 
        _incp[1]))
    _incq = gnstlib.gammainc_q(2.3, 8.3)
    print("gammainc_q(2.3, 8.3) = {0:.16e}\t err_id = {1}".format(_incq[0], 
        _incq[1]))
    _expint = gnstlib.expint(120.2, 12.2)
    print("expint(120.2, 12.2)  = {0:.16e}\t err_id = {1}".format(_expint[0], 
        _expint[1]))


def vectorized_functions():
    print("\n-------------Example vectorized functions: Gamma-------------")

    v = gnstlib.VectorDouble([1.0, 2.3, 4.3])
    r = gnstlib.VectorDouble()
    gnstlib.gamma_vec(v, r, 0)

    for i, w in enumerate(r):
        print("gamma({0}) = {1}".format(v[i], w))

    vc = gnstlib.VectorComplex([1+1j, 2+2j, 3+3j])
    rc = gnstlib.VectorComplex()
    gnstlib.gamma_vec(vc, rc, 1)

    for i, w in enumerate(rc):
        print("gamma({0}) = {1}".format(vc[i], w))


if __name__ == "__main__":
    print("GNSTLIB version 0.1 - Release October 2017")
    print("==========================================")

    basic_functions()
    gamma_functions()
    exponential_integral_functions()
    error_functions()
    incomplete_gamma_functions()
    vectorized_functions()