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

program test
    use gnstlib
    implicit none
    
    ! parameters

    integer :: err_id, i
    real(kind=dp), parameter :: x = 3.2_dp
    real(kind=dp), parameter :: xsmall = 3.2e-9_dp
    complex(kind=dp), parameter :: z = complex(3.0_dp, 4.0_dp)

    real(kind=dp) :: v(6)
    real(kind=dp), allocatable :: r(:)

    v(1:6) = (/3.2_dp, 6.3_dp, 0.4_dp, 10.5_dp, 100.6_dp, -2.7_dp/)

    err_id = 0

    write(*,*) "real functions"
    write(*,*) "x          =", x
    write(*,*) "exp(x)     =", gnstlib_exp(x)
    write(*,*) "exp2(x)    =", gnstlib_exp2(x)
    write(*,*) "exp10(x)   =", gnstlib_exp10(x)
    write(*,*) "expm1(x)   =", gnstlib_expm1(x)
    write(*,*) "expm1x(x)  =", gnstlib_expm1x(x)
    write(*,*) "expm1mx(x) =", gnstlib_expm1mx(x)
    write(*,*)
    write(*,*) "log(x)     =", gnstlib_log(x)
    write(*,*) "log2(x)    =", gnstlib_log2(x)
    write(*,*) "log10(x)   =", gnstlib_log10(x)
    write(*,*) "log1p(x)   =", gnstlib_log1p(xsmall)
    write(*,*) "log1pmx(x) =", gnstlib_log1pmx(x)
    write(*,*) "log1mexp(x)=", gnstlib_log1mexp(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "log1pexp(x)=", gnstlib_log1pexp(x)

    write(*,*) 
    write(*,*) "complex functions"
    write(*,*) "z          =", z
    write(*,*) "exp(z)     =", gnstlib_cexp(z)
    write(*,*) "log(z)     =", gnstlib_clog(z)

    write(*,*)
    write(*,*) "power functions"
    err_id = 0
    write(*,*) "powm1      =", gnstlib_powm1(-3.0_dp, -1.2_dp, err_id)
    write(*,*) "err_id     =", err_id

    write(*,*) "gamma functions"
    err_id = 0
    write(*,*) "gamma(x)      =", gnstlib_cgamma(z, err_id)
    write(*,*) "err_id        =", err_id
    write(*,*) "stirling(x)   =", gnstlib_stirling(11.5_dp, err_id)
    write(*,*) "err_id        =", err_id
    write(*,*) "gammastar(x)  =", gnstlib_gammastar(3.5_dp, err_id)
    write(*,*) "err_id        =", err_id
    write(*,*) "qgamma(x, y)  =", gnstlib_qgamma(28.5_dp, 36.2_dp, err_id)
    write(*,*) "err_id        =", err_id

    write(*,*)
    write(*,*) "vectorized gamma"
    call gnstlib_gamma_vec(6, v, r, 0)

    do i = 1, 6
        write(*,*) "gamma(",v(i),")=", r(i)
    end do

    write(*,*)
    write(*,*) "ei(x) = ", gnstlib_ei(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "e1(x) = ", gnstlib_e1(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "li(x) = ", gnstlib_li(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "ci(x) = ", gnstlib_ci(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "si(x) = ", gnstlib_si(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "chi(x) = ", gnstlib_chi(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "shi(x) = ", gnstlib_shi(x, err_id)
    write(*,*) "err_id     =", err_id

    write(*,*) "invei(x) = ", gnstlib_invei(5.4_dp, err_id)
    write(*,*) "err_id     =", err_id
    
    write(*,*)
    write(*,*) "erf(x)     =", gnstlib_erf(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "erfc(x)    =", gnstlib_erfc(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "erfcx(x)   =", gnstlib_erfcx(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "erfi(x)    =", gnstlib_erfi(x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "voigt(x)   =", gnstlib_voigt_profile(x, x, x, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "inverf(x)  =", gnstlib_inverf(0.2_dp, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "inverfc(x) =", gnstlib_inverfc(0.2_dp, err_id)
    write(*,*) "err_id     =", err_id

    write(*,*)
    write(*,*) "erf(z)     =", gnstlib_cerf(z, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "erfc(z)    =", gnstlib_cerfc(z, err_id)
    write(*,*) "err_id     =", err_id    
    write(*,*) "erfcx(z)   =", gnstlib_cerfcx(z, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "erfi(z)    =", gnstlib_cerfi(z, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "faddeeva(z)=", gnstlib_faddeeva(z, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "dawson(z)  =", gnstlib_cdawson(z, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "fresnelc(z)=", gnstlib_cfresnelc(z, err_id)
    write(*,*) "err_id     =", err_id
    write(*,*) "fresnels(z)=", gnstlib_cfresnels(z, err_id)
    write(*,*) "err_id     =", err_id

    write(*,*)
    write(*,*) "gammainc_p(a, x) = ", gnstlib_gammainc_p(2.4_dp, 1.3_dp, err_id)
    write(*,*) "gammainc_q(a, x) = ", gnstlib_gammainc_q(2.4_dp, 1.3_dp, err_id)
    write(*,*) "expint(v, x)     = ", gnstlib_expint(2.4_dp, 1.3_dp, err_id)
    write(*,*) "expint_acc(v, x) = ", gnstlib_expint_acc(2, 1e-14_dp, 1e-10_dp, err_id)

end program
