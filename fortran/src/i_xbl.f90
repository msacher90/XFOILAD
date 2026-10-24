!***********************************************************************
!   Copyright (c) 2019 D. de Vries
!   Original Copyright (c) 2000 Mark Drela
!
!   This file is part of XFoil.
!
!   XFoil is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   XFoil is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with XFoil.  If not, see <https://www.gnu.org/licenses/>.
!***********************************************************************

!*==I_XBL.f90  processed by SPAG 7.21DC at 11:25 on 11 Jan 2019
module i_xbl
    use i_blpar
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    integer, parameter :: NCOM = 73
    !
    ! Local variables
    !
    real :: amcrit, bule, cfm, cfm_d1, cfm_d2, cfm_ms, cfm_re, cfm_t1, cfm_t2, cfm_u1, cfm_u2, dwte, gambl, &
            & gm1bl, hstinv, hstinv_ms, hvrat, qinfbl, reybl, reybl_ms, reybl_re, rstbl, rstbl_ms, tkbl, &
            & tkbl_ms, xiforc, xt, xt_a1, xt_d1, xt_d2, xt_ms, xt_re, xt_t1, xt_t2, xt_u1, xt_u2, xt_x1, &
            & xt_x2, xt_xf
    real :: ampl1, ampl2, cf1, cf1_d1, cf1_ms, cf1_re, cf1_t1, cf1_u1, cf2, cf2_d2, cf2_ms, cf2_re, &
            & cf2_t2, cf2_u2, cq1, cq1_d1, cq1_ms, cq1_re, cq1_t1, cq1_u1, cq2, cq2_d2, cq2_ms, &
            & cq2_re, cq2_t2, cq2_u2, d1, d2, de1, de1_d1, de1_ms, de1_t1, de1_u1, de2, de2_d2, &
            & de2_ms, de2_t2, de2_u2, di1, di1_d1, di1_ms, di1_re, di1_s1, di1_t1, di1_u1, di2, &
            & di2_d2, di2_ms, di2_re, di2_s2, di2_t2, di2_u2, dw1, dw2, h1, h1_d1, h1_t1, h2, h2_d2, &
            & h2_t2, hc1, hc1_d1, hc1_ms, hc1_t1, hc1_u1, hc2, hc2_d2, hc2_ms, hc2_t2, hc2_u2, hk1, &
            & hk1_d1, hk1_ms, hk1_t1, hk1_u1, hk2, hk2_d2, hk2_ms, hk2_t2, hk2_u2
    real, dimension(NCOM) :: c1sav, c2sav
    real :: hs1, hs1_d1, hs1_ms, hs1_re, hs1_t1, hs1_u1, hs2, hs2_d2, hs2_ms, hs2_re, hs2_t2, &
            & hs2_u2, m1, m1_ms, m1_u1, m2, m2_ms, m2_u2, r1, r1_ms, r1_u1, r2, r2_ms, r2_u2, rt1, &
            & rt1_ms, rt1_re, rt1_t1, rt1_u1, rt2, rt2_ms, rt2_re, rt2_t2, rt2_u2, s1, s2, t1, t2, &
            & u1, u1_ms, u1_uei, u2, u2_ms, u2_uei, us1, us1_d1, us1_ms, us1_re, us1_t1, us1_u1, us2, &
            & us2_d2, us2_ms, us2_re, us2_t2, us2_u2, v1, v1_ms, v1_re, v1_u1, v2, v2_ms, v2_re, &
            & v2_u2, x1, x2
    integer :: idampv
    logical :: simi, tran, trforc, trfree, turb, wake
    real, dimension(4, 5) :: vs1, vs2
    real, dimension(4) :: vsm, vsr, vsrez, vsx
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! PARAMETER definitions
    !
    !
    ! COMMON /V_COM/
    !
    !
    ! COMMON /V_INT/
    !
    !
    ! COMMON /V_SAV/
    !
    !
    ! COMMON /V_SYS/
    !
    !
    ! COMMON /V_VAR/
    !
    !
    ! COMMON /V_VAR1/
    !
    !
    ! COMMON /V_VAR2/
    !
    !
    ! COMMON /V_VARA/
    !
contains

    subroutine blc_1from2
        implicit none
        X1 = X2
        U1 = U2
        T1 = T2
        D1 = D2
        S1 = S2
        AMPl1 = AMPl2
        U1_uei = U2_uei
        U1_ms = U2_ms
        DW1 = DW2
        H1 = H2
        H1_t1 = H2_t2
        H1_d1 = H2_d2
        M1 = M2
        M1_u1 = M2_u2
        M1_ms = M2_ms
        R1 = R2
        R1_u1 = R2_u2
        R1_ms = R2_ms
        V1 = V2
        V1_u1 = V2_u2
        V1_ms = V2_ms
        V1_re = V2_re
        HK1 = HK2
        HK1_u1 = HK2_u2
        HK1_t1 = HK2_t2
        HK1_d1 = HK2_d2
        HK1_ms = HK2_ms
        HS1 = HS2
        HS1_u1 = HS2_u2
        HS1_t1 = HS2_t2
        HS1_d1 = HS2_d2
        HS1_ms = HS2_ms
        HS1_re = HS2_re
        HC1 = HC2
        HC1_u1 = HC2_u2
        HC1_t1 = HC2_t2
        HC1_d1 = HC2_d2
        HC1_ms = HC2_ms
        RT1 = RT2
        RT1_u1 = RT2_u2
        RT1_t1 = RT2_t2
        RT1_ms = RT2_ms
        RT1_re = RT2_re
        CF1 = CF2
        CF1_u1 = CF2_u2
        CF1_t1 = CF2_t2
        CF1_d1 = CF2_d2
        CF1_ms = CF2_ms
        CF1_re = CF2_re
        DI1 = DI2
        DI1_u1 = DI2_u2
        DI1_t1 = DI2_t2
        DI1_d1 = DI2_d2
        DI1_s1 = DI2_s2
        DI1_ms = DI2_ms
        DI1_re = DI2_re
        US1 = US2
        US1_u1 = US2_u2
        US1_t1 = US2_t2
        US1_d1 = US2_d2
        US1_ms = US2_ms
        US1_re = US2_re
        CQ1 = CQ2
        CQ1_u1 = CQ2_u2
        CQ1_t1 = CQ2_t2
        CQ1_d1 = CQ2_d2
        CQ1_ms = CQ2_ms
        CQ1_re = CQ2_re
        DE1 = DE2
        DE1_u1 = DE2_u2
        DE1_t1 = DE2_t2
        DE1_d1 = DE2_d2
        DE1_ms = DE2_ms
    end subroutine blc_1from2

    subroutine blc_save1
        implicit none
        c1sav(1) = X1
        c1sav(2) = U1
        c1sav(3) = T1
        c1sav(4) = D1
        c1sav(5) = S1
        c1sav(6) = AMPl1
        c1sav(7) = U1_uei
        c1sav(8) = U1_ms
        c1sav(9) = DW1
        c1sav(10) = H1
        c1sav(11) = H1_t1
        c1sav(12) = H1_d1
        c1sav(13) = M1
        c1sav(14) = M1_u1
        c1sav(15) = M1_ms
        c1sav(16) = R1
        c1sav(17) = R1_u1
        c1sav(18) = R1_ms
        c1sav(19) = V1
        c1sav(20) = V1_u1
        c1sav(21) = V1_ms
        c1sav(22) = V1_re
        c1sav(23) = HK1
        c1sav(24) = HK1_u1
        c1sav(25) = HK1_t1
        c1sav(26) = HK1_d1
        c1sav(27) = HK1_ms
        c1sav(28) = HS1
        c1sav(29) = HS1_u1
        c1sav(30) = HS1_t1
        c1sav(31) = HS1_d1
        c1sav(32) = HS1_ms
        c1sav(33) = HS1_re
        c1sav(34) = HC1
        c1sav(35) = HC1_u1
        c1sav(36) = HC1_t1
        c1sav(37) = HC1_d1
        c1sav(38) = HC1_ms
        c1sav(39) = RT1
        c1sav(40) = RT1_u1
        c1sav(41) = RT1_t1
        c1sav(42) = RT1_ms
        c1sav(43) = RT1_re
        c1sav(44) = CF1
        c1sav(45) = CF1_u1
        c1sav(46) = CF1_t1
        c1sav(47) = CF1_d1
        c1sav(48) = CF1_ms
        c1sav(49) = CF1_re
        c1sav(50) = DI1
        c1sav(51) = DI1_u1
        c1sav(52) = DI1_t1
        c1sav(53) = DI1_d1
        c1sav(54) = DI1_s1
        c1sav(55) = DI1_ms
        c1sav(56) = DI1_re
        c1sav(57) = US1
        c1sav(58) = US1_u1
        c1sav(59) = US1_t1
        c1sav(60) = US1_d1
        c1sav(61) = US1_ms
        c1sav(62) = US1_re
        c1sav(63) = CQ1
        c1sav(64) = CQ1_u1
        c1sav(65) = CQ1_t1
        c1sav(66) = CQ1_d1
        c1sav(67) = CQ1_ms
        c1sav(68) = CQ1_re
        c1sav(69) = DE1
        c1sav(70) = DE1_u1
        c1sav(71) = DE1_t1
        c1sav(72) = DE1_d1
        c1sav(73) = DE1_ms
    end subroutine blc_save1

    subroutine blc_save2
        implicit none
        c2sav(1) = X2
        c2sav(2) = U2
        c2sav(3) = T2
        c2sav(4) = D2
        c2sav(5) = S2
        c2sav(6) = AMPl2
        c2sav(7) = U2_uei
        c2sav(8) = U2_ms
        c2sav(9) = DW2
        c2sav(10) = H2
        c2sav(11) = H2_t2
        c2sav(12) = H2_d2
        c2sav(13) = M2
        c2sav(14) = M2_u2
        c2sav(15) = M2_ms
        c2sav(16) = R2
        c2sav(17) = R2_u2
        c2sav(18) = R2_ms
        c2sav(19) = V2
        c2sav(20) = V2_u2
        c2sav(21) = V2_ms
        c2sav(22) = V2_re
        c2sav(23) = HK2
        c2sav(24) = HK2_u2
        c2sav(25) = HK2_t2
        c2sav(26) = HK2_d2
        c2sav(27) = HK2_ms
        c2sav(28) = HS2
        c2sav(29) = HS2_u2
        c2sav(30) = HS2_t2
        c2sav(31) = HS2_d2
        c2sav(32) = HS2_ms
        c2sav(33) = HS2_re
        c2sav(34) = HC2
        c2sav(35) = HC2_u2
        c2sav(36) = HC2_t2
        c2sav(37) = HC2_d2
        c2sav(38) = HC2_ms
        c2sav(39) = RT2
        c2sav(40) = RT2_u2
        c2sav(41) = RT2_t2
        c2sav(42) = RT2_ms
        c2sav(43) = RT2_re
        c2sav(44) = CF2
        c2sav(45) = CF2_u2
        c2sav(46) = CF2_t2
        c2sav(47) = CF2_d2
        c2sav(48) = CF2_ms
        c2sav(49) = CF2_re
        c2sav(50) = DI2
        c2sav(51) = DI2_u2
        c2sav(52) = DI2_t2
        c2sav(53) = DI2_d2
        c2sav(54) = DI2_s2
        c2sav(55) = DI2_ms
        c2sav(56) = DI2_re
        c2sav(57) = US2
        c2sav(58) = US2_u2
        c2sav(59) = US2_t2
        c2sav(60) = US2_d2
        c2sav(61) = US2_ms
        c2sav(62) = US2_re
        c2sav(63) = CQ2
        c2sav(64) = CQ2_u2
        c2sav(65) = CQ2_t2
        c2sav(66) = CQ2_d2
        c2sav(67) = CQ2_ms
        c2sav(68) = CQ2_re
        c2sav(69) = DE2
        c2sav(70) = DE2_u2
        c2sav(71) = DE2_t2
        c2sav(72) = DE2_d2
        c2sav(73) = DE2_ms
    end subroutine blc_save2

    subroutine blc_restore1
        implicit none
        X1 = c1sav(1)
        U1 = c1sav(2)
        T1 = c1sav(3)
        D1 = c1sav(4)
        S1 = c1sav(5)
        AMPl1 = c1sav(6)
        U1_uei = c1sav(7)
        U1_ms = c1sav(8)
        DW1 = c1sav(9)
        H1 = c1sav(10)
        H1_t1 = c1sav(11)
        H1_d1 = c1sav(12)
        M1 = c1sav(13)
        M1_u1 = c1sav(14)
        M1_ms = c1sav(15)
        R1 = c1sav(16)
        R1_u1 = c1sav(17)
        R1_ms = c1sav(18)
        V1 = c1sav(19)
        V1_u1 = c1sav(20)
        V1_ms = c1sav(21)
        V1_re = c1sav(22)
        HK1 = c1sav(23)
        HK1_u1 = c1sav(24)
        HK1_t1 = c1sav(25)
        HK1_d1 = c1sav(26)
        HK1_ms = c1sav(27)
        HS1 = c1sav(28)
        HS1_u1 = c1sav(29)
        HS1_t1 = c1sav(30)
        HS1_d1 = c1sav(31)
        HS1_ms = c1sav(32)
        HS1_re = c1sav(33)
        HC1 = c1sav(34)
        HC1_u1 = c1sav(35)
        HC1_t1 = c1sav(36)
        HC1_d1 = c1sav(37)
        HC1_ms = c1sav(38)
        RT1 = c1sav(39)
        RT1_u1 = c1sav(40)
        RT1_t1 = c1sav(41)
        RT1_ms = c1sav(42)
        RT1_re = c1sav(43)
        CF1 = c1sav(44)
        CF1_u1 = c1sav(45)
        CF1_t1 = c1sav(46)
        CF1_d1 = c1sav(47)
        CF1_ms = c1sav(48)
        CF1_re = c1sav(49)
        DI1 = c1sav(50)
        DI1_u1 = c1sav(51)
        DI1_t1 = c1sav(52)
        DI1_d1 = c1sav(53)
        DI1_s1 = c1sav(54)
        DI1_ms = c1sav(55)
        DI1_re = c1sav(56)
        US1 = c1sav(57)
        US1_u1 = c1sav(58)
        US1_t1 = c1sav(59)
        US1_d1 = c1sav(60)
        US1_ms = c1sav(61)
        US1_re = c1sav(62)
        CQ1 = c1sav(63)
        CQ1_u1 = c1sav(64)
        CQ1_t1 = c1sav(65)
        CQ1_d1 = c1sav(66)
        CQ1_ms = c1sav(67)
        CQ1_re = c1sav(68)
        DE1 = c1sav(69)
        DE1_u1 = c1sav(70)
        DE1_t1 = c1sav(71)
        DE1_d1 = c1sav(72)
        DE1_ms = c1sav(73)
    end subroutine blc_restore1

    subroutine blc_restore2
        implicit none
        X2 = c2sav(1)
        U2 = c2sav(2)
        T2 = c2sav(3)
        D2 = c2sav(4)
        S2 = c2sav(5)
        AMPl2 = c2sav(6)
        U2_uei = c2sav(7)
        U2_ms = c2sav(8)
        DW2 = c2sav(9)
        H2 = c2sav(10)
        H2_t2 = c2sav(11)
        H2_d2 = c2sav(12)
        M2 = c2sav(13)
        M2_u2 = c2sav(14)
        M2_ms = c2sav(15)
        R2 = c2sav(16)
        R2_u2 = c2sav(17)
        R2_ms = c2sav(18)
        V2 = c2sav(19)
        V2_u2 = c2sav(20)
        V2_ms = c2sav(21)
        V2_re = c2sav(22)
        HK2 = c2sav(23)
        HK2_u2 = c2sav(24)
        HK2_t2 = c2sav(25)
        HK2_d2 = c2sav(26)
        HK2_ms = c2sav(27)
        HS2 = c2sav(28)
        HS2_u2 = c2sav(29)
        HS2_t2 = c2sav(30)
        HS2_d2 = c2sav(31)
        HS2_ms = c2sav(32)
        HS2_re = c2sav(33)
        HC2 = c2sav(34)
        HC2_u2 = c2sav(35)
        HC2_t2 = c2sav(36)
        HC2_d2 = c2sav(37)
        HC2_ms = c2sav(38)
        RT2 = c2sav(39)
        RT2_u2 = c2sav(40)
        RT2_t2 = c2sav(41)
        RT2_ms = c2sav(42)
        RT2_re = c2sav(43)
        CF2 = c2sav(44)
        CF2_u2 = c2sav(45)
        CF2_t2 = c2sav(46)
        CF2_d2 = c2sav(47)
        CF2_ms = c2sav(48)
        CF2_re = c2sav(49)
        DI2 = c2sav(50)
        DI2_u2 = c2sav(51)
        DI2_t2 = c2sav(52)
        DI2_d2 = c2sav(53)
        DI2_s2 = c2sav(54)
        DI2_ms = c2sav(55)
        DI2_re = c2sav(56)
        US2 = c2sav(57)
        US2_u2 = c2sav(58)
        US2_t2 = c2sav(59)
        US2_d2 = c2sav(60)
        US2_ms = c2sav(61)
        US2_re = c2sav(62)
        CQ2 = c2sav(63)
        CQ2_u2 = c2sav(64)
        CQ2_t2 = c2sav(65)
        CQ2_d2 = c2sav(66)
        CQ2_ms = c2sav(67)
        CQ2_re = c2sav(68)
        DE2 = c2sav(69)
        DE2_u2 = c2sav(70)
        DE2_t2 = c2sav(71)
        DE2_d2 = c2sav(72)
        DE2_ms = c2sav(73)
    end subroutine blc_restore2

end module i_xbl
