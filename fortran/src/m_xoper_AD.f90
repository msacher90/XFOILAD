


module m_xoper_AD
    use s_xoper
contains

    subroutine fcpmin
        use i_xfoil
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------
        !     Finds minimum Cp on dist
        !------------------------------------------------
        !
        XCPmni = X(1)
        XCPmnv = X(1)
        CPMni = CPI(1)
        CPMnv = CPV(1)
        !
        do i = 2, N + NW
            if (CPI(i)<CPMni) then
                XCPmni = X(i)
                CPMni = CPI(i)
            endif
            if (CPV(i)<CPMnv) then
                XCPmnv = X(i)
                CPMnv = CPV(i)
            endif
        enddo
        !

        if (LVIsc) then
            CPMn = CPMnv
        else
            CPMn = CPMni
            !
            CPMnv = CPMni
            XCPmnv = XCPmni
        endif
        !
    end subroutine fcpmin




  subroutine specal
        use m_xpanel, only: qiset, ggcalc
        use i_xfoil
        use s_xfoil, only: mrcl, comset, cpcalc, clcalc, tecalc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: clm, clm1, dclm, minf_clm, msq_clm, reinf_clm
        integer :: i, irlx, itcl
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------
        !     Converges to specified alpha.
        !-----------------------------------
        !
        !---- calculate surface vorticity distributions for alpha = 0, 90 degrees
        if (.not.LGAmu .or. .not.LQAij) call ggcalc
        !
        COSa = cos(ALFa)
        SINa = sin(ALFa)
        !
        !---- superimpose suitably weighted  alpha = 0, 90  distributions
        do i = 1, N
            GAM(i) = COSa * GAMu(i, 1) + SINa * GAMu(i, 2)
            GAM_a(i) = -SINa * GAMu(i, 1) + COSa * GAMu(i, 2)
        enddo
        PSIo = COSa * GAMu(N + 1, 1) + SINa * GAMu(N + 1, 2)
        !
        call tecalc
        call qiset
        !
        !---- set initial guess for the Newton variable CLM
        clm = 1.0
        !
        !---- set corresponding  M(CLM), Re(CLM)
        call mrcl(clm, minf_clm, reinf_clm)
        call comset
        !
        !---- set corresponding CL(M)
        call clcalc(N, X, Y, GAM, GAM_a, ALFa, MINf, QINf, XCMref, YCMref, CL, CM, CDP, CL_alf, CL_msq)
        !
        !---- iterate on CLM
        do itcl = 1, 20
            !
            msq_clm = 2.0 * MINf * minf_clm
            dclm = (CL - clm) / (1.0 - CL_msq * msq_clm)
            !
            clm1 = clm
            RLX = 1.0
            !
            !------ under-relaxation loop to avoid driving M(CL) above 1
            do irlx = 1, 12
                !
                clm = clm1 + RLX * dclm
                !
                !-------- set new freestream Mach M(CLM)
                call mrcl(clm, minf_clm, reinf_clm)
                !
                !-------- if Mach is OK, go do next Newton iteration
                ! if (MATyp==1 .or. MINf==0.0 .or. minf_clm/=0.0) exit
				if (MATyp==1 .or. MINf>=0. - epsilon(0.) .and. MINf<=0. + epsilon(0.) &
				.or. .not.(minf_clm>=0. - epsilon(0.) .and. minf_clm<=0. + epsilon(0.))) exit 
                !
                RLX = 0.5 * RLX
            enddo
            !
            !------ set new CL(M)
            call comset
            call clcalc(N, X, Y, GAM, GAM_a, ALFa, MINf, QINf, XCMref, YCMref, CL, CM, CDP, CL_alf, CL_msq)
            !
            if (abs(dclm)<=1.0E-6) goto 100
            !
        enddo
        if (show_output) write (*, *) 'SPECAL:  Minf convergence failed'
        !
        !---- set final Mach, CL, Cp distributions, and hinge moment
        100  call mrcl(CL, MINf_cl, REInf_cl)
        call comset
        call clcalc(N, X, Y, GAM, GAM_a, ALFa, MINf, QINf, XCMref, YCMref, CL, CM, CDP, CL_alf, CL_msq)
        call cpcalc(N, QINv, QINf, MINf, CPI)
        if (LVIsc) then
            call cpcalc(N + NW, QVIs, QINf, MINf, CPV)
            call cpcalc(N + NW, QINv, QINf, MINf, CPI)
        else
            call cpcalc(N, QINv, QINf, MINf, CPI)
        endif
        if (LFLap) call mhinge
        !
    end subroutine specal




    subroutine speccl
        use m_xpanel, only: qiset, ggcalc
        use i_xfoil
        use s_xfoil, only: mrcl, comset, clcalc, cpcalc, tecalc
        implicit none
        real :: dalfa
        integer :: i, ital
        !-----------------------------------------
        !     Converges to specified inviscid CL.
        !-----------------------------------------
        !
        !---- calculate surface vorticity distributions for alpha = 0, 90 degrees
        if (.not.LGAmu .or. .not.LQAij) call ggcalc
        !
        !---- set freestream Mach from specified CL -- Mach will be held fixed
        call mrcl(CLSpec, MINf_cl, REInf_cl)
        call comset
        !
        !---- current alpha is the initial guess for Newton variable ALFA
        COSa = cos(ALFa)
        SINa = sin(ALFa)
        do i = 1, N
            GAM(i) = COSa * GAMu(i, 1) + SINa * GAMu(i, 2)
            GAM_a(i) = -SINa * GAMu(i, 1) + COSa * GAMu(i, 2)
        enddo
        PSIo = COSa * GAMu(N + 1, 1) + SINa * GAMu(N + 1, 2)
        !
        !---- get corresponding CL, CL_alpha, CL_Mach
        call clcalc(N, X, Y, GAM, GAM_a, ALFa, MINf, QINf, XCMref, YCMref, CL, CM, CDP, CL_alf, CL_msq)
        !
        !---- Newton loop for alpha to get specified inviscid CL
        do ital = 1, 20
            !
            dalfa = (CLSpec - CL) / CL_alf
            RLX = 1.0
            !
            ALFa = ALFa + RLX * dalfa
            !
            !------ set new surface speed distribution
            COSa = cos(ALFa)
            SINa = sin(ALFa)
            do i = 1, N
                GAM(i) = COSa * GAMu(i, 1) + SINa * GAMu(i, 2)
                GAM_a(i) = -SINa * GAMu(i, 1) + COSa * GAMu(i, 2)
            enddo
            PSIo = COSa * GAMu(N + 1, 1) + SINa * GAMu(N + 1, 2)
            !
            !------ set new CL(alpha)
            call clcalc(N, X, Y, GAM, GAM_a, ALFa, MINf, QINf, XCMref, YCMref, CL, CM, CDP, CL_alf, CL_msq)
            !
            if (abs(dalfa)<=1.0E-6) goto 100
        enddo
        if (show_output) write (*, *) 'SPECCL:  CL convergence failed'
        !
        !---- set final surface speed and Cp distributions
        100  call tecalc
        call qiset
        if (LVIsc) then
            call cpcalc(N + NW, QVIs, QINf, MINf, CPV)
            call cpcalc(N + NW, QINv, QINf, MINf, CPI)
        else
            call cpcalc(N, QINv, QINf, MINf, CPI)
        endif
        !
    end subroutine speccl



    function viscal(Niter1)
        use m_xpanel, only: qiset, qvfue, iblpan, qwcalc, uicalc, xicalc, gamqv, stfind, stmove, xywake, qdcalc
        use m_xbl, only: update, setbl
        use m_xsolve, only: blsolv
        use m_userio, only: aski
        use i_xfoil
        use s_xbl, only: iblsys
        use s_xfoil, only: mrcl, comset, cdcalc, cpcalc, clcalc
        implicit none

        logical :: viscal
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: Niter1
        ! intent (in) Niter1
        !
        ! Local variables
        !
        real :: cdpdif, ddef, delp, dpds, duds, edef, fnum, hfrac, hk, hki, hkm, hkmax, patt, pdef, pdefi, &
                & pdefm, psep
        integer :: i, ibl, iblm, iblp, idec, ione, is, iten, iter, izero, lu, niter
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        !
        ! Local variables
        !
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------
        !     Converges viscous operating point
        !----------------------------------------
        !
        !---- convergence tolerance
        viscal = .true.
        !
        niter = Niter1
        !
        !---- calculate wake trajectory from current inviscid solution if necessary
        if (.not.LWAke) call xywake
        !
        !---- set velocities on wake from airfoil vorticity for alpha=0, 90
        call qwcalc
        !
        !---- set velocities on airfoil and wake for initial alpha
        call qiset
        !
        if (.not.LIPan) then
            !
            if (LBLini) call gamqv
            !
            !----- locate stagnation point arc length position and panel index
            call stfind
            !
            !----- set  BL position -> panel position  pointers
            call iblpan
            !
            !----- calculate surface arc length array for current stagnation point location
            call xicalc
            !
            !----- set  BL position -> system line  pointers
            call iblsys
            !
        endif
        !
        !---- set inviscid BL edge velocity UINV from QINV
        call uicalc
        !
        if (.not.LBLini) then
            !
            !----- set initial Ue from inviscid Ue
            do ibl = 1, NBL(1)
                UEDg(ibl, 1) = UINv(ibl, 1)
            enddo
            !
            do ibl = 1, NBL(2)
                UEDg(ibl, 2) = UINv(ibl, 2)
            enddo
            !
        endif
        !
        if (LVConv) then
            !----- set correct CL if converged point exists
            call qvfue
            if (LVIsc) then
                call cpcalc(N + NW, QVIs, QINf, MINf, CPV)
                call cpcalc(N + NW, QINv, QINf, MINf, CPI)
            else
                call cpcalc(N, QINv, QINf, MINf, CPI)
            endif
            call gamqv
            call clcalc(N, X, Y, GAM, GAM_a, ALFa, MINf, QINf, XCMref, YCMref, CL, CM, CDP, CL_alf, CL_msq)
            call cdcalc
        endif
        !
        !---- set up source influence matrix if it doesn't exist
        if (.not.LWDij .or. .not.LADij) call qdcalc
        !
        !---- Newton iteration for entire BL solution
        if (niter==0) call aski('Enter number of iterations^', niter)
        if (show_output) then
            write (*, *)
            write (*, *) 'Solving BL system ...'
        endif
        do iter = 1, niter
            !
            !------ fill Newton system for BL variables
            viscal = setbl()
            if (abort_on_nan .and. .not. viscal) return
            !
            !------ solve Newton system with custom solver
            call blsolv
            !
            !------ update BL variables
            call update
            !
            if (LALfa) then
                !------- set new freestream Mach, Re from new CL
                call mrcl(CL, MINf_cl, REInf_cl)
                call comset
            else
                !------- set new inviscid speeds QINV and UINV for new alpha
                call qiset
                call uicalc
            endif
            !
            !------ calculate edge velocities QVIS(.) from UEDG(..)
            call qvfue
            !
            !------ set GAM distribution from QVIS
            call gamqv
            !
            !------ relocate stagnation point
            call stmove
            !
            !------ set updated CL,CD
            call clcalc(N, X, Y, GAM, GAM_a, ALFa, MINf, QINf, XCMref, YCMref, CL, CM, CDP, CL_alf, CL_msq)
            call cdcalc
            !
            !------ display changes and test for convergence
            if (RLX<1.0 .and. show_output) write (*, 99001) iter, RMSbl, RMXbl, VMXbl, IMXbl, ISMxbl, RLX
            !....................................................................
            99001 format (/1x, i3, '   rms: ', e10.4, '   max: ', e10.4, 3x, a1, ' at ', i4, i3, '   RLX:', f6.3)
            ! if (RLX==1.0 .and. show_output) write (*, 99002) iter, RMSbl, RMXbl, VMXbl, IMXbl, ISMxbl
			if (RLX>=1. - epsilon(1.) .and. RLX<=1. + epsilon(1.) &
			.and. show_output) write (*, 99002) iter, RMSbl, RMXbl, VMXbl, IMXbl, ISMxbl 
            99002 format (/1x, i3, '   rms: ', e10.4, '   max: ', e10.4, 3x, a1, ' at ', i4, i3)
            cdpdif = CD - CDF
            if (show_output) write (*, 99003) ALFa / DTOr, CL, CM, CD, CDF, cdpdif
            99003 format (1x, 3x, '   a =', f7.3, '      CL =', f8.4/1x, 3x, '  Cm =', f8.4, &
                    '     CD =', f9.5, '   =>   CDf =', f9.5, &
                    &'    CDp =', f9.5)
            !         CDSURF = CDP + CDF
            !         WRITE(*,2025) CDSURF, CDF, CDP

            if (RMSbl<xad_eps1 .and. iter>=xad_minit) then
                LVConv = .true.
                AVIsc = ALFa
                MVIsc = MINf
                goto 100
            endif
            !
        enddo
        if (show_output) write (*, *) 'VISCAL:  Convergence failed'
        !
        100  call cpcalc(N + NW, QINv, QINf, MINf, CPI)
        call cpcalc(N + NW, QVIs, QINf, MINf, CPV)
        if (LFLap) call mhinge

        ! is = 1
        ! hkmax = 0.
        ! hkm = 0.0
        ! psep = 0.
        ! patt = 0.
        ! do ibl = 2, IBLte(is)
            ! hki = DSTr(ibl, is) / THEt(ibl, is)
            ! hkmax = max(hki, hkmax)
            ! if (hkm<4.0 .and. hki>=4.0) then
                ! hfrac = (4.0 - hkm) / (hki - hkm)
                ! pdefm = UEDg(ibl - 1, is)**2 * THEt(ibl - 1, is)
                ! pdefi = UEDg(ibl, is)**2 * THEt(ibl, is)
                ! psep = pdefm * (1.0 - hfrac) + pdefi * hfrac
            ! endif
            ! if (hkm>4.0 .and. hki<4.0) then
                ! hfrac = (4.0 - hkm) / (hki - hkm)
                ! pdefm = UEDg(ibl - 1, is)**2 * THEt(ibl - 1, is)
                ! pdefi = UEDg(ibl, is)**2 * THEt(ibl, is)
                ! patt = pdefm * (1.0 - hfrac) + pdefi * hfrac
            ! endif
            ! hkm = hki
        ! enddo
        ! delp = patt - psep

        ! if (show_output) write (*, 99004) ACRit(is), hkmax, CD, 2.0 * psep, 2.0 * patt, 2.0 * delp, XOCtr(is)
        ! 99004 format (1x, f10.3, f10.4, f11.6, 3F11.6, f10.4, '     #')

        ! izero = ichar('0')

        !! c      fnum = acrit(is)
        ! fnum = XSTrip(is) * 100.0

        ! iten = int(fnum / 9.99999)
        ! ione = int((fnum - float(10 * iten)) / 0.99999)
        ! idec = int((fnum - float(10 * iten) - float(ione)) / 0.09999)

        ! FNAme = char(iten + izero) // char(ione + izero) // char(idec + izero) // '.bl'
        ! lu = 44
        ! open (lu, file = FNAme, status = 'unknown')
        ! rewind (lu)
        ! write (lu, '(a,a)') '#       s         ue          H          P         K ', '        x    -m du/dx'
              !! 1234567890 1234567890 1234567890 1234567890 1234567890 1234567890
        ! do ibl = 2, IBLte(is)
            ! iblm = max(ibl - 1, 2)
            ! iblp = min(ibl + 1, IBLte(is))
            ! i = IPAn(ibl, is)
            ! hk = DSTr(ibl, is) / THEt(ibl, is)
            ! ddef = DSTr(ibl, is) * UEDg(ibl, is)
            ! pdef = THEt(ibl, is) * UEDg(ibl, is)**2
            ! edef = TSTr(ibl, is) * UEDg(ibl, is)**3 * 0.5
            ! duds = (UEDg(iblp, is) - UEDg(iblm, is)) / (XSSi(iblp, is) - XSSi(iblm, is))
            ! dpds = -ddef * duds
            ! write (lu, 99005) XSSi(ibl, is), UEDg(ibl, is), hk, pdef, edef, X(i), dpds
            ! 99005 format (1x, 3F11.4, 2F11.6, f11.3, e14.6)
        ! enddo
        ! close (lu)

    end function viscal


end module m_xoper_AD