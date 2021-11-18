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

%module gnstlib

%{
  /* the resulting C file should be built as a python extension */
  #define SWIG_FILE_WITH_INIT
  /* include modules */
  #include "../gnstlib/include/gnstlib_constants.hpp"
  #include "../gnstlib/include/gnstlib_basic.hpp"
  #include "../gnstlib/include/gnstlib.hpp" 
%}

/******************************* INCLUDES *************************************/
%include <std_complex.i>
%include <std_vector.i>

namespace std {
  %template (VectorDouble) vector<double>;
  %template (VectorComplex) vector<std::complex<double>>;
};

/******************************* TYPES ****************************************/
// return list [value, err_id]
%apply int &OUTPUT { int & err_id };
%apply std::complex<double> &OUTPUT { std::complex<double> & si };
%apply std::complex<double> &OUTPUT { std::complex<double> & ci };


// GNSTLIB machine precision and mathematical constants
%include "../gnstlib/include/gnstlib_constants.hpp"
// GNSTLIB basic functions
%include "../gnstlib/include/gnstlib_basic.hpp"
// GNSTLIB special functions
%include "../gnstlib/include/gnstlib.hpp"
