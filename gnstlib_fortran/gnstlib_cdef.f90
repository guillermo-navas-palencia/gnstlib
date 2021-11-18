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

! C functions declaration
interface
    function gnstlib_exp_c(x) bind(C, name="gnstlib_exp")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_exp_c
        real(c_double), intent(in), value :: x
    end function gnstlib_exp_c

    function gnstlib_exp2_c(x) bind(C, name="gnstlib_exp2")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_exp2_c
        real(c_double), intent(in), value :: x
    end function gnstlib_exp2_c

    function gnstlib_exp10_c(x) bind(C, name="gnstlib_exp10")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_exp10_c
        real(c_double), intent(in), value :: x
    end function gnstlib_exp10_c

    function gnstlib_expm1_c(x) bind(C, name="gnstlib_expm1")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_expm1_c
        real(c_double), intent(in), value :: x
    end function gnstlib_expm1_c

    function gnstlib_expm1x_c(x) bind(C, name="gnstlib_expm1x")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_expm1x_c
        real(c_double), intent(in), value :: x
    end function gnstlib_expm1x_c

    function gnstlib_expm1mx_c(x) bind(C, name="gnstlib_expm1mx")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_expm1mx_c
        real(c_double), intent(in), value :: x
    end function gnstlib_expm1mx_c

    function gnstlib_log_c(x) bind(C, name="gnstlib_log")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_log_c
        real(c_double), intent(in), value :: x
    end function gnstlib_log_c

    function gnstlib_log2_c(x) bind(C, name="gnstlib_log2")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_log2_c
        real(c_double), intent(in), value :: x
    end function gnstlib_log2_c

    function gnstlib_log10_c(x) bind(C, name="gnstlib_log10")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_log10_c
        real(c_double), intent(in), value :: x
    end function gnstlib_log10_c

    function gnstlib_log1p_c(x) bind(C, name="gnstlib_log1p")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_log1p_c
        real(c_double), intent(in), value :: x
    end function gnstlib_log1p_c

    function gnstlib_log1pmx_c(x) bind(C, name="gnstlib_log1pmx")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_log1pmx_c
        real(c_double), intent(in), value :: x
    end function gnstlib_log1pmx_c

    function gnstlib_log1mexp_c(x, err_id) bind(C, name="gnstlib_log1mexp")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_log1mexp_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_log1mexp_c

    function gnstlib_log1pexp_c(x) bind(C, name="gnstlib_log1pexp")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_log1pexp_c
        real(c_double), intent(in), value :: x
    end function gnstlib_log1pexp_c

    function gnstlib_cexp_c(z) bind(C, name="gnstlib_cexp")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cexp_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_cexp_c

    function gnstlib_clog_c(z) bind(C, name="gnstlib_clog")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_clog_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_clog_c

    function gnstlib_cos_c(x) bind(C, name="gnstlib_cos")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_cos_c
        real(c_double), intent(in), value :: x
    end function gnstlib_cos_c

    function gnstlib_sin_c(x) bind(C, name="gnstlib_sin")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_sin_c
        real(c_double), intent(in), value :: x
    end function gnstlib_sin_c

    function gnstlib_tan_c(x) bind(C, name="gnstlib_tan")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_tan_c
        real(c_double), intent(in), value :: x
    end function gnstlib_tan_c

    function gnstlib_sec_c(x) bind(C, name="gnstlib_sec")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_sec_c
        real(c_double), intent(in), value :: x
    end function gnstlib_sec_c

    function gnstlib_csc_c(x) bind(C, name="gnstlib_csc")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_csc_c
        real(c_double), intent(in), value :: x
    end function gnstlib_csc_c

    function gnstlib_cot_c(x) bind(C, name="gnstlib_cot")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_cot_c
        real(c_double), intent(in), value :: x
    end function gnstlib_cot_c

    function gnstlib_acos_c(x) bind(C, name="gnstlib_acos")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_acos_c
        real(c_double), intent(in), value :: x
    end function gnstlib_acos_c

    function gnstlib_asin_c(x) bind(C, name="gnstlib_asin")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_asin_c
        real(c_double), intent(in), value :: x
    end function gnstlib_asin_c

    function gnstlib_atan_c(x) bind(C, name="gnstlib_atan")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_atan_c
        real(c_double), intent(in), value :: x
    end function gnstlib_atan_c

    function gnstlib_atan2_c(y, x) bind(C, name="gnstlib_atan2")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_atan2_c
        real(c_double), intent(in), value :: y, x
    end function gnstlib_atan2_c

    function gnstlib_asec_c(x) bind(C, name="gnstlib_asec")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_asec_c
        real(c_double), intent(in), value :: x
    end function gnstlib_asec_c

    function gnstlib_acsc_c(x) bind(C, name="gnstlib_acsc")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_acsc_c
        real(c_double), intent(in), value :: x
    end function gnstlib_acsc_c

    function gnstlib_acot_c(x) bind(C, name="gnstlib_acot")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_acot_c
        real(c_double), intent(in), value :: x
    end function gnstlib_acot_c

    function gnstlib_cospi_c(x) bind(C, name="gnstlib_cospi")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_cospi_c
        real(c_double), intent(in), value :: x
    end function gnstlib_cospi_c

    function gnstlib_sinpi_c(x) bind(C, name="gnstlib_sinpi")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_sinpi_c
        real(c_double), intent(in), value :: x
    end function gnstlib_sinpi_c

    function gnstlib_ccos_c(z) bind(C, name="gnstlib_ccos")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_ccos_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_ccos_c

    function gnstlib_csin_c(z) bind(C, name="gnstlib_csin")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_csin_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_csin_c

    function gnstlib_ctan_c(z) bind(C, name="gnstlib_ctan")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_ctan_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_ctan_c

    function gnstlib_csec_c(z) bind(C, name="gnstlib_csec")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_csec_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_csec_c

    function gnstlib_ccsc_c(z) bind(C, name="gnstlib_ccsc")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_ccsc_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_ccsc_c

    function gnstlib_ccot_c(z) bind(C, name="gnstlib_ccot")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_ccot_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_ccot_c

    function gnstlib_cacos_c(z) bind(C, name="gnstlib_cacos")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cacos_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_cacos_c

    function gnstlib_casin_c(z) bind(C, name="gnstlib_casin")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_casin_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_casin_c

    function gnstlib_catan_c(z) bind(C, name="gnstlib_catan")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_catan_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_catan_c

    function gnstlib_casec_c(z) bind(C, name="gnstlib_casec")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_casec_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_casec_c

    function gnstlib_cacsc_c(z) bind(C, name="gnstlib_cacsc")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cacsc_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_cacsc_c

    function gnstlib_cacot_c(z) bind(C, name="gnstlib_cacot")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cacot_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_cacot_c

    function gnstlib_ccospi_c(z) bind(C, name="gnstlib_ccospi")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_ccospi_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_ccospi_c

    function gnstlib_csinpi_c(z) bind(C, name="gnstlib_csinpi")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_csinpi_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_csinpi_c

    function gnstlib_cosh_c(x) bind(C, name="gnstlib_cosh")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_cosh_c
        real(c_double), intent(in), value :: x
    end function gnstlib_cosh_c

    function gnstlib_sinh_c(x) bind(C, name="gnstlib_sinh")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_sinh_c
        real(c_double), intent(in), value :: x
    end function gnstlib_sinh_c

    function gnstlib_tanh_c(x) bind(C, name="gnstlib_tanh")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_tanh_c
        real(c_double), intent(in), value :: x
    end function gnstlib_tanh_c

    function gnstlib_sech_c(x) bind(C, name="gnstlib_sech")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_sech_c
        real(c_double), intent(in), value :: x
    end function gnstlib_sech_c

    function gnstlib_csch_c(x) bind(C, name="gnstlib_csch")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_csch_c
        real(c_double), intent(in), value :: x
    end function gnstlib_csch_c

    function gnstlib_coth_c(x) bind(C, name="gnstlib_coth")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_coth_c
        real(c_double), intent(in), value :: x
    end function gnstlib_coth_c

    function gnstlib_acosh_c(x) bind(C, name="gnstlib_acosh")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_acosh_c
        real(c_double), intent(in), value :: x
    end function gnstlib_acosh_c

    function gnstlib_asinh_c(x) bind(C, name="gnstlib_asinh")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_asinh_c
        real(c_double), intent(in), value :: x
    end function gnstlib_asinh_c

    function gnstlib_atanh_c(x) bind(C, name="gnstlib_atanh")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_atanh_c
        real(c_double), intent(in), value :: x
    end function gnstlib_atanh_c

    function gnstlib_asech_c(x) bind(C, name="gnstlib_asech")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_asech_c
        real(c_double), intent(in), value :: x
    end function gnstlib_asech_c

    function gnstlib_acsch_c(x) bind(C, name="gnstlib_acsch")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_acsch_c
        real(c_double), intent(in), value :: x
    end function gnstlib_acsch_c

    function gnstlib_acoth_c(x) bind(C, name="gnstlib_acoth")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_acoth_c
        real(c_double), intent(in), value :: x
    end function gnstlib_acoth_c

    function gnstlib_ccosh_c(z) bind(C, name="gnstlib_ccosh")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_ccosh_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_ccosh_c

    function gnstlib_csinh_c(z) bind(C, name="gnstlib_csinh")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_csinh_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_csinh_c

    function gnstlib_ctanh_c(z) bind(C, name="gnstlib_ctanh")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_ctanh_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_ctanh_c

    function gnstlib_csech_c(z) bind(C, name="gnstlib_csech")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_csech_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_csech_c

    function gnstlib_ccsch_c(z) bind(C, name="gnstlib_ccsch")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_ccsch_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_ccsch_c

    function gnstlib_ccoth_c(z) bind(C, name="gnstlib_ccoth")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_ccoth_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_ccoth_c

    function gnstlib_cacosh_c(z) bind(C, name="gnstlib_cacosh")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cacosh_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_cacosh_c

    function gnstlib_casinh_c(z) bind(C, name="gnstlib_casinh")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_casinh_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_casinh_c

    function gnstlib_catanh_c(z) bind(C, name="gnstlib_catanh")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_catanh_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_catanh_c

    function gnstlib_casech_c(z) bind(C, name="gnstlib_casech")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_casech_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_casech_c

    function gnstlib_cacsch_c(z) bind(C, name="gnstlib_cacsch")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cacsch_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_cacsch_c

    function gnstlib_cacoth_c(z) bind(C, name="gnstlib_cacoth")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cacoth_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_cacoth_c

    function gnstlib_pow_c(x, y) bind(C, name="gnstlib_pow")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_pow_c
        real(c_double), intent(in), value :: x, y
    end function gnstlib_pow_c

    function gnstlib_sqrt_c(x) bind(C, name="gnstlib_sqrt")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_sqrt_c
        real(c_double), intent(in), value :: x
    end function gnstlib_sqrt_c

    function gnstlib_hypot_c(x, y) bind(C, name="gnstlib_hypot")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_hypot_c
        real(c_double), intent(in), value :: x, y
    end function gnstlib_hypot_c

    function gnstlib_cbrt_c(x) bind(C, name="gnstlib_cbrt")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_cbrt_c
        real(c_double), intent(in), value :: x
    end function gnstlib_cbrt_c

    function gnstlib_powm1_c(x, y, err_id) bind(C, name="gnstlib_powm1")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_powm1_c
        real(c_double), intent(in), value :: x, y
        integer(c_int), intent(out) :: err_id
    end function gnstlib_powm1_c

    function gnstlib_cpow_c(x, y) bind(C, name="gnstlib_cpow")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cpow_c
        complex(c_double_complex), intent(in), value :: x, y
    end function gnstlib_cpow_c

    function gnstlib_csqrt_c(z) bind(C, name="gnstlib_csqrt")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_csqrt_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_csqrt_c
    
    function gnstlib_gamma_c(x, err_id) bind(C, name="gnstlib_gamma")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_gamma_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_gamma_c

    function gnstlib_gammaln_c(x, err_id) bind(C, name="gnstlib_gammaln")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_gammaln_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_gammaln_c

    function gnstlib_gammasign_c(x) bind(C, name="gnstlib_gammasign")
        use iso_c_binding
        implicit none
        integer(c_int) :: gnstlib_gammasign_c
        real(c_double), intent(in), value :: x
    end function gnstlib_gammasign_c

    function gnstlib_factorial_c(n, err_id) bind(C, name="gnstlib_factorial")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_factorial_c
        integer(c_int), intent(in), value :: n
        integer(c_int), intent(out) :: err_id
    end function gnstlib_factorial_c

    function gnstlib_qgamma_c(x, y, err_id) bind(C, name="gnstlib_qgamma")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_qgamma_c
        real(c_double), intent(in), value :: x, y
        integer(c_int), intent(out) :: err_id
    end function gnstlib_qgamma_c

    function gnstlib_cgamma_c(z, err_id) bind(C, name="gnstlib_cgamma")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cgamma_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cgamma_c

    function gnstlib_loggamma_c(z, err_id) bind(C, name="gnstlib_loggamma")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_loggamma_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_loggamma_c

    function gnstlib_rgamma_c(z) bind(C, name="gnstlib_rgamma")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_rgamma_c
        complex(c_double_complex), intent(in), value :: z
    end function gnstlib_rgamma_c

    function gnstlib_stirling_c(x, err_id) bind(C, name="gnstlib_stirling")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_stirling_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_stirling_c

    function gnstlib_gammastar_c(x, err_id) bind(C, name="gnstlib_gammastar")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_gammastar_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_gammastar_c

    function gnstlib_auxgam_c(x, err_id) bind(C, name="gnstlib_auxgam")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_auxgam_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_auxgam_c

    subroutine gnstlib_gamma_vec_c(n, v, r, option) bind(C, name="gnstlib_gamma_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_gamma_vec_c

    subroutine gnstlib_gammaln_vec_c(n, v, r, option) bind(C, name="gnstlib_gammaln_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_gammaln_vec_c

    subroutine gnstlib_cgamma_vec_c(n, v, r, option) bind(C, name="gnstlib_cgamma_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cgamma_vec_c

    subroutine gnstlib_loggamma_vec_c(n, v, r, option) bind(C, name="gnstlib_loggamma_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_loggamma_vec_c

    subroutine gnstlib_rgamma_vec_c(n, v, r, option) bind(C, name="gnstlib_rgamma_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_rgamma_vec_c

    function gnstlib_ei_c(x, err_id) bind(C, name="gnstlib_ei")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_ei_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_ei_c

    function gnstlib_e1_c(x, err_id) bind(C, name="gnstlib_e1")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_e1_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_e1_c

    function gnstlib_li_c(x, err_id) bind(C, name="gnstlib_li")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_li_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_li_c

    function gnstlib_ci_c(x, err_id) bind(C, name="gnstlib_ci")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_ci_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_ci_c

    function gnstlib_si_c(x, err_id) bind(C, name="gnstlib_si")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_si_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_si_c

    function gnstlib_chi_c(x, err_id) bind(C, name="gnstlib_chi")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_chi_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_chi_c

    function gnstlib_shi_c(x, err_id) bind(C, name="gnstlib_shi")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_shi_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_shi_c

    function gnstlib_cei_c(z, err_id) bind(C, name="gnstlib_cei")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cei_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cei_c

    function gnstlib_ce1_c(z, err_id) bind(C, name="gnstlib_ce1")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_ce1_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_ce1_c

    function gnstlib_cli_c(z, err_id) bind(C, name="gnstlib_cli")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cli_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cli_c

    function gnstlib_cci_c(z, err_id) bind(C, name="gnstlib_cci")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cci_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cci_c

    function gnstlib_csi_c(z, err_id) bind(C, name="gnstlib_csi")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_csi_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_csi_c

    function gnstlib_cchi_c(z, err_id) bind(C, name="gnstlib_cchi")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cchi_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cchi_c

    function gnstlib_cshi_c(z, err_id) bind(C, name="gnstlib_cshi")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cshi_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cshi_c

    function gnstlib_invei_c(x, err_id) bind(C, name="gnstlib_invei")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_invei_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_invei_c

    function gnstlib_invli_c(x, err_id) bind(C, name="gnstlib_invli")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_invli_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_invli_c

    function gnstlib_erf_c(x, err_id) bind(C, name="gnstlib_erf")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_erf_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_erf_c

    subroutine gnstlib_ei_vec_c(n, v, r, option) bind(C, name="gnstlib_ei_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_ei_vec_c

    subroutine gnstlib_e1_vec_c(n, v, r, option) bind(C, name="gnstlib_e1_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_e1_vec_c

    subroutine gnstlib_li_vec_c(n, v, r, option) bind(C, name="gnstlib_li_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_li_vec_c

    subroutine gnstlib_ci_vec_c(n, v, r, option) bind(C, name="gnstlib_ci_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_ci_vec_c

    subroutine gnstlib_si_vec_c(n, v, r, option) bind(C, name="gnstlib_si_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_si_vec_c

    subroutine gnstlib_chi_vec_c(n, v, r, option) bind(C, name="gnstlib_chi_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_chi_vec_c

    subroutine gnstlib_shi_vec_c(n, v, r, option) bind(C, name="gnstlib_shi_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_shi_vec_c

    subroutine gnstlib_cei_vec_c(n, v, r, option) bind(C, name="gnstlib_cei_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cei_vec_c

    subroutine gnstlib_ce1_vec_c(n, v, r, option) bind(C, name="gnstlib_ce1_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_ce1_vec_c

    subroutine gnstlib_cli_vec_c(n, v, r, option) bind(C, name="gnstlib_cli_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cli_vec_c

    subroutine gnstlib_cci_vec_c(n, v, r, option) bind(C, name="gnstlib_cci_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cci_vec_c

    subroutine gnstlib_csi_vec_c(n, v, r, option) bind(C, name="gnstlib_csi_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_csi_vec_c

    subroutine gnstlib_cchi_vec_c(n, v, r, option) bind(C, name="gnstlib_cchi_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cchi_vec_c

    subroutine gnstlib_cshi_vec_c(n, v, r, option) bind(C, name="gnstlib_cshi_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cshi_vec_c

    function gnstlib_erfc_c(x, err_id) bind(C, name="gnstlib_erfc")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_erfc_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_erfc_c

    function gnstlib_erfcx_c(x, err_id) bind(C, name="gnstlib_erfcx")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_erfcx_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_erfcx_c

    function gnstlib_erfi_c(x, err_id) bind(C, name="gnstlib_erfi")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_erfi_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_erfi_c

    function gnstlib_dawson_c(x, err_id) bind(C, name="gnstlib_dawson")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_dawson_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_dawson_c

    function gnstlib_fresnelc_c(x, err_id) bind(C, name="gnstlib_fresnelc")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_fresnelc_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_fresnelc_c

    function gnstlib_fresnels_c(x, err_id) bind(C, name="gnstlib_fresnels")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_fresnels_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_fresnels_c

    function gnstlib_voigt_profile_c(x, sigma, gamma, err_id) bind(C, name="gnstlib_voigt_profile")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_voigt_profile_c
        real(c_double), intent(in), value :: x, sigma, gamma
        integer(c_int), intent(out) :: err_id
    end function gnstlib_voigt_profile_c

    function gnstlib_cerf_c(z, err_id) bind(C, name="gnstlib_cerf")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cerf_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cerf_c

    function gnstlib_cerfc_c(z, err_id) bind(C, name="gnstlib_cerfc")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cerfc_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cerfc_c

    function gnstlib_cerfcx_c(z, err_id) bind(C, name="gnstlib_cerfcx")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cerfcx_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cerfcx_c

    function gnstlib_cerfi_c(z, err_id) bind(C, name="gnstlib_cerfi")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cerfi_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cerfi_c

    function gnstlib_cdawson_c(z, err_id) bind(C, name="gnstlib_cdawson")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cdawson_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cdawson_c

    function gnstlib_faddeeva_c(z, err_id) bind(C, name="gnstlib_faddeeva")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_faddeeva_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_faddeeva_c

    function gnstlib_cfresnelc_c(z, err_id) bind(C, name="gnstlib_cfresnelc")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cfresnelc_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cfresnelc_c

    function gnstlib_cfresnels_c(z, err_id) bind(C, name="gnstlib_cfresnels")
        use iso_c_binding
        implicit none
        complex(c_double_complex) :: gnstlib_cfresnels_c
        complex(c_double_complex), intent(in), value :: z
        integer(c_int), intent(out) :: err_id
    end function gnstlib_cfresnels_c

    function gnstlib_inverf_c(x, err_id) bind(C, name="gnstlib_inverf")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_inverf_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_inverf_c

    function gnstlib_inverfc_c(x, err_id) bind(C, name="gnstlib_inverfc")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_inverfc_c
        real(c_double), intent(in), value :: x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_inverfc_c

    subroutine gnstlib_erf_vec_c(n, v, r, option) bind(C, name="gnstlib_erf_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_erf_vec_c

    subroutine gnstlib_erfc_vec_c(n, v, r, option) bind(C, name="gnstlib_erfc_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_erfc_vec_c

    subroutine gnstlib_erfcx_vec_c(n, v, r, option) bind(C, name="gnstlib_erfcx_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_erfcx_vec_c

    subroutine gnstlib_erfi_vec_c(n, v, r, option) bind(C, name="gnstlib_erfi_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_erfi_vec_c

    subroutine gnstlib_dawson_vec_c(n, v, r, option) bind(C, name="gnstlib_dawson_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_dawson_vec_c

    subroutine gnstlib_fresnelc_vec_c(n, v, r, option) bind(C, name="gnstlib_fresnelc_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_fresnelc_vec_c

    subroutine gnstlib_fresnels_vec_c(n, v, r, option) bind(C, name="gnstlib_fresnels_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_fresnels_vec_c

    subroutine gnstlib_inverf_vec_c(n, v, r, option) bind(C, name="gnstlib_inverf_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_inverf_vec_c

    subroutine gnstlib_inverfc_vec_c(n, v, r, option) bind(C, name="gnstlib_inverfc_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        real(c_double), intent(in)  :: v(n)
        real(c_double), intent(out) :: r(n)
    end subroutine gnstlib_inverfc_vec_c

    subroutine gnstlib_cerf_vec_c(n, v, r, option) bind(C, name="gnstlib_cerf_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cerf_vec_c

    subroutine gnstlib_cerfc_vec_c(n, v, r, option) bind(C, name="gnstlib_cerfc_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cerfc_vec_c

    subroutine gnstlib_cerfcx_vec_c(n, v, r, option) bind(C, name="gnstlib_cerfcx_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cerfcx_vec_c

    subroutine gnstlib_cerfi_vec_c(n, v, r, option) bind(C, name="gnstlib_cerfi_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cerfi_vec_c

    subroutine gnstlib_cdawson_vec_c(n, v, r, option) bind(C, name="gnstlib_cdawson_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cdawson_vec_c

    subroutine gnstlib_faddeeva_vec_c(n, v, r, option) bind(C, name="gnstlib_faddeeva_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_faddeeva_vec_c

    subroutine gnstlib_cfresnelc_vec_c(n, v, r, option) bind(C, name="gnstlib_cfresnelc_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cfresnelc_vec_c

    subroutine gnstlib_cfresnels_vec_c(n, v, r, option) bind(C, name="gnstlib_cfresnels_vec")
        use iso_c_binding
        implicit none

        integer(c_int), intent(in), value :: n, option
        complex(c_double_complex), intent(in)  :: v(n)
        complex(c_double_complex), intent(out) :: r(n)
    end subroutine gnstlib_cfresnels_vec_c

    function gnstlib_gammainc_p_c(a, x, err_id) bind(C, name="gnstlib_gammainc_p")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_gammainc_p_c
        real(c_double), intent(in), value :: a, x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_gammainc_p_c

    function gnstlib_gammainc_q_c(a, x, err_id) bind(C, name="gnstlib_gammainc_q")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_gammainc_q_c
        real(c_double), intent(in), value :: a, x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_gammainc_q_c

    function gnstlib_expint_c(v, x, err_id) bind(C, name="gnstlib_expint")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_expint_c
        real(c_double), intent(in), value :: v, x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_expint_c

    function gnstlib_expint_acc_c(n, e, x, err_id) bind(C, name="gnstlib_expint_acc")
        use iso_c_binding
        implicit none
        real(c_double) :: gnstlib_expint_acc_c
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in), value :: e, x
        integer(c_int), intent(out) :: err_id
    end function gnstlib_expint_acc_c

end interface