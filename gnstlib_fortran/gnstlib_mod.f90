! Copyright (C) 2017 A. Gil, G. Navas-Palencia, J. Segura and N. M. Temme

! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:

! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
! OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
! WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

module gnstlib

    ! include interface declaration for external C API
    include "gnstlib_cdef.f90"

    ! double precision
    integer, parameter :: dp = kind(1.d0)

contains

    ! *********************************************************************
    ! Exponential and logarithmic functions
    ! *********************************************************************     

    function gnstlib_exp(x)
    !   Compute exp(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_exp
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_exp = gnstlib_exp_c(x)
    end function gnstlib_exp

    function gnstlib_exp2(x)
    !   Compute exp2(x) = 2**x for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_exp2
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_exp2 = gnstlib_exp2_c(x)
    end function gnstlib_exp2

    function gnstlib_exp10(x)
    !   Compute exp10(x) = 10**x for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_exp10
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_exp10 = gnstlib_exp10_c(x)
    end function gnstlib_exp10

    function gnstlib_expm1(x)
    !   Compute exp(x) - 1 for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_expm1
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_expm1 = gnstlib_expm1_c(x)
    end function gnstlib_expm1

    function gnstlib_expm1x(x)
    !   Compute (exp(x) - 1) / x for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_expm1x
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_expm1x = gnstlib_expm1x_c(x)
    end function gnstlib_expm1x

    function gnstlib_expm1mx(x)
    !   Compute (exp(x) - 1 - x) / (0.5 * x**2) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_expm1mx
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_expm1mx = gnstlib_expm1mx_c(x)
    end function gnstlib_expm1mx

    function gnstlib_log(x)
    !   Compute log(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_log
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_log = gnstlib_log_c(x)
    end function gnstlib_log

    function gnstlib_log2(x)
    !   Compute log2(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_log2
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_log2 = gnstlib_log2_c(x)
    end function gnstlib_log2

    function gnstlib_log10(x)
    !   Compute log10(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_log10
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_log10 = gnstlib_log10_c(x)
    end function gnstlib_log10

    function gnstlib_log1p(x)
    !   Compute log1p(x) = log(1 + x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_log1p
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_log1p = gnstlib_log1p_c(x)
    end function gnstlib_log1p

    function gnstlib_log1pmx(x)
    !   Compute log1pmx(x) = log(1 + x) - x for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_log1pmx
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_log1pmx = gnstlib_log1pmx_c(x)
    end function gnstlib_log1pmx

    function gnstlib_log1pexp(x)
    !   Compute log1pexp(x) = log(1 + exp(x)) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_log1pexp
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_log1pexp = gnstlib_log1pexp_c(x)
    end function gnstlib_log1pexp

    function gnstlib_log1mexp(x, err_id)
    !   Compute log1mexp(x) = log(1 - exp(x)) for real x
    !   Constraint: x < 0.0
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_log1mexp
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_log1mexp = gnstlib_log1mexp_c(x, err_id)
    end function gnstlib_log1mexp

    function gnstlib_clog(z)
    !   Compute log(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_clog
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_clog = gnstlib_clog_c(z)
    end function gnstlib_clog

    function gnstlib_cexp(z)
    !   Compute exp(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cexp
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_cexp = gnstlib_cexp_c(z)
    end function gnstlib_cexp

    ! *********************************************************************
    ! Trigonometric functions
    ! *********************************************************************

    function gnstlib_cos(x)
    !   Compute cos(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_cos
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_cos = gnstlib_cos_c(x)
    end function gnstlib_cos

    function gnstlib_sin(x)
    !   Compute sin(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_sin
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_sin = gnstlib_sin_c(x)
    end function gnstlib_sin

    function gnstlib_tan(x)
    !   Compute tan(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_tan
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_tan = gnstlib_tan_c(x)
    end function gnstlib_tan

    function gnstlib_sec(x)
    !   Compute sec(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_sec
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_sec = gnstlib_sec_c(x)
    end function gnstlib_sec

    function gnstlib_csc(x)
    !   Compute csc(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_csc
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_csc = gnstlib_csc_c(x)
    end function gnstlib_csc

    function gnstlib_cot(x)
    !   Compute cot(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_cot
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_cot = gnstlib_cot_c(x)
    end function gnstlib_cot

    function gnstlib_acos(x)
    !   Compute acos(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_acos
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_acos = gnstlib_acos_c(x)
    end function gnstlib_acos

    function gnstlib_asin(x)
    !   Compute asin(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_asin
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_asin = gnstlib_asin_c(x)
    end function gnstlib_asin

    function gnstlib_atan(x)
    !   Compute atan(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_atan
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_atan = gnstlib_atan_c(x)
    end function gnstlib_atan

    function gnstlib_atan2(y, x)
    !   Compute atan2(x) for real y and x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_atan2
    !   scalar argument
        real(kind=dp), intent(in) :: y, x
        
        gnstlib_atan2 = gnstlib_atan2_c(y, x)
    end function gnstlib_atan2

    function gnstlib_asec(x)
    !   Compute asec(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_asec
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_asec = gnstlib_asec_c(x)
    end function gnstlib_asec

    function gnstlib_acsc(x)
    !   Compute acsc(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_acsc
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_acsc = gnstlib_acsc_c(x)
    end function gnstlib_acsc

    function gnstlib_acot(x)
    !   Compute acot(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_acot
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_acot = gnstlib_acot_c(x)
    end function gnstlib_acot

    function gnstlib_cospi(x)
    !   Compute cos(pi * x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_cospi
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_cospi = gnstlib_cospi_c(x)
    end function gnstlib_cospi

    function gnstlib_sinpi(x)
    !   Compute sin(pi * x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_sinpi
    !   scalar argument
        real(kind=dp), intent(in) :: x
        
        gnstlib_sinpi = gnstlib_sinpi_c(x)
    end function gnstlib_sinpi

    function gnstlib_ccos(z)
    !   Compute cos(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_ccos
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_ccos = gnstlib_ccos_c(z)
    end function gnstlib_ccos

    function gnstlib_csin(z)
    !   Compute sin(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_csin
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_csin = gnstlib_csin_c(z)
    end function gnstlib_csin

    function gnstlib_ctan(z)
    !   Compute tan(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_ctan
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_ctan = gnstlib_ctan_c(z)
    end function gnstlib_ctan

    function gnstlib_csec(z)
    !   Compute sec(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_csec
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_csec = gnstlib_csec_c(z)
    end function gnstlib_csec

    function gnstlib_ccsc(z)
    !   Compute csc(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_ccsc
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_ccsc = gnstlib_ccsc_c(z)
    end function gnstlib_ccsc

    function gnstlib_ccot(z)
    !   Compute cot(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_ccot
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_ccot = gnstlib_ccot_c(z)
    end function gnstlib_ccot

    function gnstlib_cacos(z)
    !   Compute acos(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cacos
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_cacos = gnstlib_cacos_c(z)
    end function gnstlib_cacos

    function gnstlib_casin(z)
    !   Compute asin(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_casin
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_casin = gnstlib_casin_c(z)
    end function gnstlib_casin

    function gnstlib_catan(z)
    !   Compute atan(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_catan
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_catan = gnstlib_catan_c(z)
    end function gnstlib_catan

    function gnstlib_casec(z)
    !   Compute asec(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_casec
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_casec = gnstlib_casec_c(z)
    end function gnstlib_casec

    function gnstlib_cacsc(z)
    !   Compute acsc(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cacsc
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_cacsc = gnstlib_cacsc_c(z)
    end function gnstlib_cacsc

    function gnstlib_cacot(z)
    !   Compute acot(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cacot
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_cacot = gnstlib_cacot_c(z)
    end function gnstlib_cacot

    function gnstlib_ccospi(z)
    !   Compute cos(pi * z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_ccospi
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_ccospi = gnstlib_ccospi_c(z)
    end function gnstlib_ccospi

    function gnstlib_csinpi(z)
    !   Compute sin(pi * z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_csinpi
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_csinpi = gnstlib_csinpi_c(z)
    end function gnstlib_csinpi

    ! *********************************************************************
    ! Hyperbolic functions
    ! *********************************************************************

    function gnstlib_cosh(x)
    !   Compute cosh(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_cosh
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_cosh = gnstlib_cosh_c(x)
    end function gnstlib_cosh

    function gnstlib_sinh(x)
    !   Compute sinh(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_sinh
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_sinh = gnstlib_sinh_c(x)
    end function gnstlib_sinh

    function gnstlib_tanh(x)
    !   Compute tanh(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_tanh
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_tanh = gnstlib_tanh_c(x)
    end function gnstlib_tanh

    function gnstlib_sech(x)
    !   Compute sech(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_sech
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_sech = gnstlib_sech_c(x)
    end function gnstlib_sech

    function gnstlib_csch(x)
    !   Compute sech(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_csch
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_csch = gnstlib_csch_c(x)
    end function gnstlib_csch

    function gnstlib_coth(x)
    !   Compute coth(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_coth
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_coth = gnstlib_coth_c(x)
    end function gnstlib_coth

    function gnstlib_acosh(x)
    !   Compute acosh(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_acosh
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_acosh = gnstlib_acosh_c(x)
    end function gnstlib_acosh

    function gnstlib_asinh(x)
    !   Compute asinh(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_asinh
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_asinh = gnstlib_asinh_c(x)
    end function gnstlib_asinh

    function gnstlib_atanh(x)
    !   Compute atanh(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_atanh
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_atanh = gnstlib_atanh_c(x)
    end function gnstlib_atanh

    function gnstlib_asech(x)
    !   Compute asech(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_asech
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_asech = gnstlib_asech_c(x)
    end function gnstlib_asech

    function gnstlib_acsch(x)
    !   Compute acsch(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_acsch
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_acsch = gnstlib_acsch_c(x)
    end function gnstlib_acsch

    function gnstlib_acoth(x)
    !   Compute acoth(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_acoth
    !   scalar argument
        real(kind=dp), intent(in) :: x

        gnstlib_acoth = gnstlib_acoth_c(x)
    end function gnstlib_acoth

    function gnstlib_ccosh(z)
    !   Compute cosh(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_ccosh
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_ccosh = gnstlib_ccosh_c(z)
    end function gnstlib_ccosh

    function gnstlib_csinh(z)
    !   Compute sinh(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_csinh
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_csinh = gnstlib_csinh_c(z)
    end function gnstlib_csinh

    function gnstlib_ctanh(z)
    !   Compute tanh(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_ctanh
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_ctanh = gnstlib_ctanh_c(z)
    end function gnstlib_ctanh

    function gnstlib_csech(z)
    !   Compute sech(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_csech
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_csech = gnstlib_csech_c(z)
    end function gnstlib_csech

    function gnstlib_ccsch(z)
    !   Compute csch(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_ccsch
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_ccsch = gnstlib_ccsch_c(z)
    end function gnstlib_ccsch

    function gnstlib_ccoth(z)
    !   Compute coth(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_ccoth
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_ccoth = gnstlib_ccoth_c(z)
    end function gnstlib_ccoth

    function gnstlib_cacosh(z)
    !   Compute acosh(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cacosh
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_cacosh = gnstlib_cacosh_c(z)
    end function gnstlib_cacosh

    function gnstlib_casinh(z)
    !   Compute asinh(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_casinh
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_casinh = gnstlib_casinh_c(z)
    end function gnstlib_casinh

    function gnstlib_catanh(z)
    !   Compute atanh(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_catanh
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_catanh = gnstlib_catanh_c(z)
    end function gnstlib_catanh

    function gnstlib_casech(z)
    !   Compute asech(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_casech
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_casech = gnstlib_casech_c(z)
    end function gnstlib_casech

    function gnstlib_cacsch(z)
    !   Compute acsch(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cacsch
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_cacsch = gnstlib_cacsch_c(z)
    end function gnstlib_cacsch

    function gnstlib_cacoth(z)
    !   Compute acoth(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cacoth
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_cacoth = gnstlib_cacoth_c(z)
    end function gnstlib_cacoth

    ! *********************************************************************
    ! Power functions
    ! *********************************************************************

    function gnstlib_pow(x, y)
    !   Compute x^y for real x and y

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_pow
    !   scalar arguments
        real(kind=dp), intent(in) :: x, y

        gnstlib_pow = gnstlib_pow_c(x, y)
    end function gnstlib_pow

    function gnstlib_sqrt(x)
    !   Compute sqrt(x) for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_sqrt
    !   scalar arguments
        real(kind=dp), intent(in) :: x

        gnstlib_sqrt = gnstlib_sqrt_c(x)
    end function gnstlib_sqrt

    function gnstlib_hypot(x, y)
    !   Compute hypot(x, y) = sqrt(x^2 + y^2)

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_hypot
    !   scalar arguments
        real(kind=dp), intent(in) :: x, y

        gnstlib_hypot = gnstlib_hypot_c(x, y)
    end function gnstlib_hypot

    function gnstlib_cbrt(x)
    !   Compute cbrt(x) = cubic root for real x

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_cbrt
    !   scalar arguments
        real(kind=dp), intent(in) :: x

        gnstlib_cbrt = gnstlib_cbrt_c(x)
    end function gnstlib_cbrt

    function gnstlib_powm1(x, y, err_id)
    !   Compute x^y - 1 for real x and y

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_powm1
    !   scalar arguments
        real(kind=dp), intent(in) :: x, y
        integer, intent(out)      :: err_id

        gnstlib_powm1 = gnstlib_powm1_c(x, y, err_id)
    end function gnstlib_powm1

    function gnstlib_cpow(x, y)
    !   Compute x^y for complex z and/or y

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cpow
    !   scalar arguments
        complex(kind=dp), intent(in) :: x, y

        gnstlib_cpow = gnstlib_cpow_c(x, y)
    end function gnstlib_cpow

    function gnstlib_csqrt(z)
    !   Compute sqrt(z) for complex z or x < 0.0

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_csqrt
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_csqrt = gnstlib_csqrt_c(z)
    end function gnstlib_csqrt

    ! *********************************************************************
    ! Gamma functions
    ! *********************************************************************

    function gnstlib_gamma(x, err_id)
    !   Compute Gamma(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_gamma
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_gamma = gnstlib_gamma_c(x, err_id)
    end function gnstlib_gamma

    function gnstlib_gammaln(x, err_id)
    !   Compute Gammaln(x) = Log(Gamma(|x|)) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_gammaln
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_gammaln = gnstlib_gammaln_c(x, err_id)
    end function gnstlib_gammaln

    function gnstlib_gammasign(x)
    !   Compute sign of Gamma(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_gammasign
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        
        gnstlib_gammasign = gnstlib_gammasign_c(x)
    end function gnstlib_gammasign

    function gnstlib_factorial(n, err_id)
    !   Compute n! for integer n, n > 0
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_factorial
    !   scalar arguments
        integer, intent(in)  :: n
        integer, intent(out) :: err_id
        
        gnstlib_factorial = gnstlib_factorial_c(n, err_id)
    end function gnstlib_factorial

    function gnstlib_qgamma(x, y, err_id)
    !   Compute quotient of Gamma(x)/Gamma(y) for real x and y
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_qgamma
    !   scalar arguments
        real(kind=dp), intent(in) :: x, y
        integer, intent(out)      :: err_id
        
        gnstlib_qgamma = gnstlib_qgamma_c(x, y, err_id)
    end function gnstlib_qgamma

    function gnstlib_cgamma(z, err_id)
    !   Compute Gamma(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cgamma
    !   scalar argument
        complex(kind=dp), intent(in) :: z
        integer, intent(out)      :: err_id

        gnstlib_cgamma = gnstlib_cgamma_c(z, err_id)
    end function gnstlib_cgamma

    function gnstlib_loggamma(z, err_id)
    !   Compute LogGamma(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_loggamma
    !   scalar argument
        complex(kind=dp), intent(in) :: z
        integer, intent(out)      :: err_id

        gnstlib_loggamma = gnstlib_loggamma_c(z, err_id)
    end function gnstlib_loggamma

    function gnstlib_rgamma(z)
    !   Compute 1/ Gamma(z) for complex z

    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_rgamma
    !   scalar argument
        complex(kind=dp), intent(in) :: z

        gnstlib_rgamma = gnstlib_rgamma_c(z)
    end function gnstlib_rgamma

    function gnstlib_stirling(x, err_id)
    !   Compute Stirling series S(x) for real x > 0
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_stirling
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_stirling = gnstlib_stirling_c(x, err_id)
    end function gnstlib_stirling

    function gnstlib_gammastar(x, err_id)
    !   Compute regulated function Gamma*(x) for real x > 0
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_gammastar
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_gammastar = gnstlib_gammastar_c(x, err_id)
    end function gnstlib_gammastar

    function gnstlib_auxgam(x, err_id)
    !   Compute auxiliary function g(x) in 1/Gamma(1+x) = 1 + x(x-1)g(x) for
    !   x in [-1, 1]

    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_auxgam
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id

        gnstlib_auxgam = gnstlib_auxgam_c(x, err_id)
    end function gnstlib_auxgam

    subroutine gnstlib_gamma_vec(n, v, r, option)
    !   Compute Gamma(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_gamma_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_gamma_vec

    subroutine gnstlib_gammaln_vec(n, v, r, option)
    !   Compute Gammaln(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_gammaln_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_gammaln_vec

    subroutine gnstlib_cgamma_vec(n, v, r, option)
    !   Compute Gamma(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cgamma_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cgamma_vec

    subroutine gnstlib_loggamma_vec(n, v, r, option)
    !   Compute LogGamma(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_loggamma_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_loggamma_vec

    subroutine gnstlib_rgamma_vec(n, v, r, option)
    !   Compute 1/Gamma(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_rgamma_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_rgamma_vec

    ! *********************************************************************
    ! Exponential, Logarithmic, trigonometric and hyperbolic integrals
    ! *********************************************************************

    function gnstlib_ei(x, err_id)
    !   Compute Ei(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_ei
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_ei = gnstlib_ei_c(x, err_id)
    end function gnstlib_ei

    function gnstlib_e1(x, err_id)
    !   Compute E1(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_e1
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_e1 = gnstlib_e1_c(x, err_id)
    end function gnstlib_e1

    function gnstlib_li(x, err_id)
    !   Compute li(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_li
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_li = gnstlib_li_c(x, err_id)
    end function gnstlib_li

    function gnstlib_ci(x, err_id)
    !   Compute Ci(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_ci
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_ci = gnstlib_ci_c(x, err_id)
    end function gnstlib_ci

    function gnstlib_si(x, err_id)
    !   Compute Si(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_si
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_si = gnstlib_si_c(x, err_id)
    end function gnstlib_si

    function gnstlib_chi(x, err_id)
    !   Compute Chi(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_chi
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_chi = gnstlib_chi_c(x, err_id)
    end function gnstlib_chi

    function gnstlib_shi(x, err_id)
    !   Compute Shi(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_shi
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_shi = gnstlib_shi_c(x, err_id)
    end function gnstlib_shi

    function gnstlib_cei(z, err_id)
    !   Compute Ei(x) for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cei
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cei = gnstlib_cei_c(z, err_id)
    end function gnstlib_cei

    function gnstlib_ce1(z, err_id)
    !   Compute E1(x) for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_ce1
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_ce1 = gnstlib_ce1_c(z, err_id)
    end function gnstlib_ce1

    function gnstlib_cli(z, err_id)
    !   Compute li(x) for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cli
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cli = gnstlib_cli_c(z, err_id)
    end function gnstlib_cli

    function gnstlib_cci(z, err_id)
    !   Compute Ci(x) for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cci
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cci = gnstlib_cci_c(z, err_id)
    end function gnstlib_cci

    function gnstlib_csi(z, err_id)
    !   Compute Si(x) for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_csi
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_csi = gnstlib_csi_c(z, err_id)
    end function gnstlib_csi

    function gnstlib_cchi(z, err_id)
    !   Compute Chi(x) for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cchi
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cchi = gnstlib_cchi_c(z, err_id)
    end function gnstlib_cchi

    function gnstlib_cshi(z, err_id)
    !   Compute Shi(x) for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cshi
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cshi = gnstlib_cshi_c(z, err_id)
    end function gnstlib_cshi

    function gnstlib_invei(x, err_id)
    !   Compute inverse of Ei(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_invei
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_invei = gnstlib_invei_c(x, err_id)
    end function gnstlib_invei

    function gnstlib_invli(x, err_id)
    !   Compute inverse of li(x) for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_invli
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_invli = gnstlib_invli_c(x, err_id)
    end function gnstlib_invli

    subroutine gnstlib_ei_vec(n, v, r, option)
    !   Compute Ei(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_ei_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_ei_vec

    subroutine gnstlib_e1_vec(n, v, r, option)
    !   Compute E1(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_e1_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_e1_vec

    subroutine gnstlib_li_vec(n, v, r, option)
    !   Compute li(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_li_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_li_vec

    subroutine gnstlib_ci_vec(n, v, r, option)
    !   Compute Ci(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_ci_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_ci_vec

    subroutine gnstlib_si_vec(n, v, r, option)
    !   Compute Si(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_si_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_si_vec

    subroutine gnstlib_chi_vec(n, v, r, option)
    !   Compute Chi(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_chi_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_chi_vec

    subroutine gnstlib_shi_vec(n, v, r, option)
    !   Compute Shi(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_shi_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_shi_vec

    subroutine gnstlib_cei_vec(n, v, r, option)
    !   Compute Ei(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cei_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cei_vec

    subroutine gnstlib_ce1_vec(n, v, r, option)
    !   Compute E1(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_ce1_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_ce1_vec

    subroutine gnstlib_cli_vec(n, v, r, option)
    !   Compute li(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cli_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cli_vec

    subroutine gnstlib_cci_vec(n, v, r, option)
    !   Compute Ci(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cci_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cci_vec

    subroutine gnstlib_csi_vec(n, v, r, option)
    !   Compute Si(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_csi_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_csi_vec

    subroutine gnstlib_cchi_vec(n, v, r, option)
    !   Compute Chi(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cchi_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cchi_vec

    subroutine gnstlib_cshi_vec(n, v, r, option)
    !   Compute Shi(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cshi_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cshi_vec

    ! *********************************************************************
    ! Error functions, Dawson's and Fresnel Integrals
    ! *********************************************************************

    function gnstlib_erf(x, err_id)
    !   Compute error function for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_erf
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_erf = gnstlib_erf_c(x, err_id)
    end function gnstlib_erf

    function gnstlib_erfc(x, err_id)
    !   Compute complementary error function for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_erfc
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_erfc = gnstlib_erfc_c(x, err_id)
    end function gnstlib_erfc

    function gnstlib_erfcx(x, err_id)
    !   Compute scaled complementary error function for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_erfcx
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_erfcx = gnstlib_erfcx_c(x, err_id)
    end function gnstlib_erfcx

    function gnstlib_erfi(x, err_id)
    !   Compute imaginary error function for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_erfi
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_erfi = gnstlib_erfi_c(x, err_id)
    end function gnstlib_erfi

    function gnstlib_dawson(x, err_id)
    !   Compute Dawson's integral for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_dawson
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_dawson = gnstlib_dawson_c(x, err_id)
    end function gnstlib_dawson

    function gnstlib_fresnelc(x, err_id)
    !   Compute Fresnel C integral for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_fresnelc
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_fresnelc = gnstlib_fresnelc_c(x, err_id)
    end function gnstlib_fresnelc

    function gnstlib_fresnels(x, err_id)
    !   Compute Fresnel S integral for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_fresnels
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_fresnels = gnstlib_fresnels_c(x, err_id)
    end function gnstlib_fresnels

    function gnstlib_voigt_profile(x, sigma, gamma, err_id)
    !   Compute Voigt profile
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_voigt_profile
    !   scalar arguments
        real(kind=dp), intent(in) :: x, sigma, gamma
        integer, intent(out)      :: err_id
        
        gnstlib_voigt_profile = gnstlib_voigt_profile_c(x, sigma, gamma, err_id)
    end function gnstlib_voigt_profile

    function gnstlib_cerf(z, err_id)
    !   Compute error function for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cerf
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cerf = gnstlib_cerf_c(z, err_id)
    end function gnstlib_cerf

    function gnstlib_cerfc(z, err_id)
    !   Compute complementary error function for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cerfc
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cerfc = gnstlib_cerfc_c(z, err_id)
    end function gnstlib_cerfc

    function gnstlib_cerfcx(z, err_id)
    !   Compute scaled complemtary error function for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cerfcx
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cerfcx = gnstlib_cerfcx_c(z, err_id)
    end function gnstlib_cerfcx

    function gnstlib_cerfi(z, err_id)
    !   Compute imaginary error function for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cerfi
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cerfi = gnstlib_cerfi_c(z, err_id)
    end function gnstlib_cerfi

    function gnstlib_cdawson(z, err_id)
    !   Compute Dawson's integral for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cdawson
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cdawson = gnstlib_cdawson_c(z, err_id)
    end function gnstlib_cdawson

    function gnstlib_faddeeva(z, err_id)
    !   Compute Faddeeva for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_faddeeva
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_faddeeva = gnstlib_faddeeva_c(z, err_id)
    end function gnstlib_faddeeva

    function gnstlib_cfresnelc(z, err_id)
    !   Compute Fresnel C integral for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cfresnelc
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cfresnelc = gnstlib_cfresnelc_c(z, err_id)
    end function gnstlib_cfresnelc

    function gnstlib_cfresnels(z, err_id)
    !   Compute Fresnel S integral for complex z
    
    !   implicit none statement
        implicit none
    !   function return value
        complex(kind=dp)             :: gnstlib_cfresnels
    !   scalar arguments
        complex(kind=dp), intent(in) :: z
        integer, intent(out)         :: err_id
        
        gnstlib_cfresnels = gnstlib_cfresnels_c(z, err_id)
    end function gnstlib_cfresnels

    function gnstlib_inverf(x, err_id)
    !   Compute inverse error function for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_inverf
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_inverf = gnstlib_inverf_c(x, err_id)
    end function gnstlib_inverf

    function gnstlib_inverfc(x, err_id)
    !   Compute inverse complementary error function for real x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_inverfc
    !   scalar arguments
        real(kind=dp), intent(in) :: x
        integer, intent(out)      :: err_id
        
        gnstlib_inverfc = gnstlib_inverfc_c(x, err_id)
    end function gnstlib_inverfc

    subroutine gnstlib_erf_vec(n, v, r, option)
    !   Compute erf(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_erf_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_erf_vec

    subroutine gnstlib_erfc_vec(n, v, r, option)
    !   Compute erfc(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_erfc_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_erfc_vec

    subroutine gnstlib_erfcx_vec(n, v, r, option)
    !   Compute erfcx(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_erfcx_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_erfcx_vec

    subroutine gnstlib_erfi_vec(n, v, r, option)
    !   Compute erfi(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_erfi_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_erfi_vec

    subroutine gnstlib_dawson_vec(n, v, r, option)
    !   Compute dawson(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_dawson_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_dawson_vec

    subroutine gnstlib_fresnelc_vec(n, v, r, option)
    !   Compute fresnelc(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_fresnelc_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_fresnelc_vec

    subroutine gnstlib_fresnels_vec(n, v, r, option)
    !   Compute fresnels(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_fresnels_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_fresnels_vec

    subroutine gnstlib_inverf_vec(n, v, r, option)
    !   Compute inverf(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_inverf_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_inverf_vec

    subroutine gnstlib_inverfc_vec(n, v, r, option)
    !   Compute inverfc(array) for real arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                     :: n, option
    !   Array arguments
        real(kind=dp), intent(in)               :: v(n)
        real(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                 :: ierr
    !   Intrinsic
        intrinsic                               :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_inverfc_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_inverfc_vec

    subroutine gnstlib_cerf_vec(n, v, r, option)
    !   Compute erf(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cerf_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cerf_vec

    subroutine gnstlib_cerfc_vec(n, v, r, option)
    !   Compute erfc(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cerfc_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cerfc_vec

    subroutine gnstlib_cerfcx_vec(n, v, r, option)
    !   Compute erfcx(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cerfcx_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cerfcx_vec

    subroutine gnstlib_cerfi_vec(n, v, r, option)
    !   Compute erfi(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cerfi_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cerfi_vec

    subroutine gnstlib_cdawson_vec(n, v, r, option)
    !   Compute dawson(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cdawson_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cdawson_vec

    subroutine gnstlib_faddeeva_vec(n, v, r, option)
    !   Compute faddeeva(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_faddeeva_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_faddeeva_vec

    subroutine gnstlib_cfresnelc_vec(n, v, r, option)
    !   Compute fresnelc(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cfresnelc_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cfresnelc_vec

    subroutine gnstlib_cfresnels_vec(n, v, r, option)
    !   Compute fresnels(array) for complex arrays
    
    !   implicit none statement
        implicit none
    !   Scalar arguments
        integer, intent(in)                        :: n, option
    !   Array arguments
        complex(kind=dp), intent(in)               :: v(n)
        complex(kind=dp), allocatable, intent(out) :: r(:)
    !   Local scalar
        integer                                    :: ierr
    !   Intrinsic
        intrinsic                                  :: allocated

        ierr = 0

        ! clean initial allocation if allocated
        if (allocated(r)) then 
            deallocate(r)
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        else
            allocate(r(n), stat=ierr)
            if (ierr /= 0) then
                go to 100
            end if
        end if

        call gnstlib_cfresnels_vec_c(n, v, r, option)

100     continue
        return
    end subroutine gnstlib_cfresnels_vec
    
    ! *********************************************************************
    ! Incomplete gamma and generalized exponential integral
    ! *********************************************************************

    function gnstlib_gammainc_p(a, x, err_id)
    !   Compute regularized incomplete gamma function P(a, x) for real a and x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_gammainc_p
    !   scalar arguments
        real(kind=dp), intent(in) :: a, x
        integer, intent(out)      :: err_id
        
        gnstlib_gammainc_p = gnstlib_gammainc_p_c(a, x, err_id)
    end function gnstlib_gammainc_p

    function gnstlib_gammainc_q(a, x, err_id)
    !   Compute regularized incomplete gamma function Q(a, x) for real a and x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_gammainc_q
    !   scalar arguments
        real(kind=dp), intent(in) :: a, x
        integer, intent(out)      :: err_id
        
        gnstlib_gammainc_q = gnstlib_gammainc_q_c(a, x, err_id)
    end function gnstlib_gammainc_q

    function gnstlib_expint(v, x, err_id)
    !   Compute generalized exponential integral Ev(x) for real v and x
    
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_expint
    !   scalar arguments
        real(kind=dp), intent(in) :: v, x
        integer, intent(out)      :: err_id
        
        gnstlib_expint = gnstlib_expint_c(v, x, err_id)
    end function gnstlib_expint

    function gnstlib_expint_acc(n, e, x, err_id)
    !   Compute accurate generalized exponential integral Ev(x) for real v and x
        
    !   implicit none statement
        implicit none
    !   function return value
        real(kind=dp)             :: gnstlib_expint_acc
    !   scalar arguments
        integer, intent(in)       :: n
        real(kind=dp), intent(in) :: e, x
        integer, intent(out)      :: err_id
        
        gnstlib_expint_acc = gnstlib_expint_acc_c(n, e, x, err_id)
    end function gnstlib_expint_acc

end module gnstlib