!***********************************************************************
!   XFOILAD - C interface of the adjoint-enabled XFOIL library
!
!   All reals are C doubles, all arrays are contiguous (Fortran order).
!   Logical flags are passed as C ints (0 / 1).
!***********************************************************************
module xad_api
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    implicit none
    !> panel nodes of the last primal analysis, valid until an adjoint pass
    !> (xad_keepgeom: the next analysis on the same nodes keeps the geometry)
    real(c_double), allocatable, save :: geo_x(:), geo_y(:)
    logical, save :: geo_ok = .false.
contains


    !> Resets XFOIL to its default state (must be called once after loading).
    subroutine api_init() bind(c, name='xad_init')
        use xad_core_diff, only: xad_init
        call xad_init()
        geo_ok = .false.
    end subroutine api_init


    !> Flow conditions. re <= 0 gives an inviscid analysis.
    subroutine api_set_flow(re, mach, ncrit_top, ncrit_bot, xtr_top, xtr_bot, verbose) &
            bind(c, name='xad_set_flow')
        use i_xfoil_diff
        real(c_double), intent(in) :: re, mach, ncrit_top, ncrit_bot, xtr_top, xtr_bot
        integer(c_int), intent(in) :: verbose

        REInf1 = re
        LVIsc = re > 0.0
        MINf1 = mach
        ACRit(1) = ncrit_top
        ACRit(2) = ncrit_bot
        XSTrip(1) = xtr_top
        XSTrip(2) = xtr_bot
        show_output = verbose /= 0
        LBLini = .false.
        LVConv = .false.
    end subroutine api_set_flow


    !> Viscous Newton convergence tolerance (XFOIL default 1e-4).
    subroutine api_set_tolerance(eps1) bind(c, name='xad_set_tolerance')
        use i_xfoil_diff, only: xad_eps1
        real(c_double), intent(in) :: eps1
        xad_eps1 = eps1
    end subroutine api_set_tolerance


    !> Minimum distance of the stagnation point from the panel nodes, as a
    !> fraction of the panel length (0.01 by default).
    subroutine api_set_stag_frac(f) bind(c, name='xad_set_stag_frac')
        use i_xfoil_diff, only: STFrac
        real(c_double), intent(in) :: f
        STFrac = f
    end subroutine api_set_stag_frac


    !> Initial alpha (deg) used as Newton guess in specified-CL mode.
    subroutine api_set_alpha_guess(adeg) bind(c, name='xad_set_alpha_guess')
        use i_xfoil_diff, only: ALFa, DTOr
        real(c_double), intent(in) :: adeg
        ALFa = adeg * DTOr
    end subroutine api_set_alpha_guess


    !> Primal analysis of one operating point.
    !>   imode = 1 : op = alpha (deg)      imode = 2 : op = CL
    !>   warm  = 0 : boundary layers restarted from scratch
    !>   warm  = 1 : restart from the current BL state (continuation)
    !>   res(1:7) = cl, cd, cdf, cm, alpha, xtr_top, xtr_bot
    !>   cp(n), cf(n) = nodal pressure and skin friction coefficients
    subroutine api_analyze(n, x, y, op, imode, niter, warm, res, cp, cf, conv) &
            bind(c, name='xad_analyze')
        use xad_core_diff, only: xad_eval
        use i_xfoil_diff, only: XOCtr, xad_minit, REInf1, ACRit, xad_keepgeom
        integer(c_int), intent(in) :: n, imode, niter, warm
        real(c_double), intent(in) :: x(n), y(n), op
        real(c_double), intent(out) :: res(7), cp(n), cf(n)
        integer(c_int), intent(out) :: conv
        real :: cl, cd, cdf, cm, alfa, rew, ncw(2)
        logical :: lconv

        xad_minit = 1
        rew = REInf1
        ncw = ACRit
        xad_keepgeom = .false.
        if (geo_ok) then
            if (size(geo_x) == n) xad_keepgeom = all(geo_x == x) .and. all(geo_y == y)
        endif
        call xad_eval(n, x, y, op, rew, ncw, imode, niter, warm /= 0, cl, cd, cdf, cm, alfa, cp, cf, lconv)
        xad_keepgeom = .false.
        if (allocated(geo_x)) deallocate(geo_x, geo_y)
        allocate(geo_x(n), geo_y(n))
        geo_x = x
        geo_y = y
        geo_ok = .true.
        res(1) = cl
        res(2) = cd
        res(3) = cdf
        res(4) = cm
        res(5) = alfa
        res(6) = XOCtr(1)
        res(7) = XOCtr(2)
        conv = merge(1, 0, lconv)
    end subroutine api_analyze


    !> Adjoint (reverse) sweep around the converged state of the last
    !> api_analyze call made with the same (x, y, op, imode).
    !>   seed(1:5)  = weights of cl, cd, cdf, cm, alpha in the scalar output
    !>   cp_seed(n), cf_seed(n) = weights of the nodal Cp and Cf
    !>   gx, gy, gop = gradient of  sum(seed*outputs) + sum(cp_seed*cp) + sum(cf_seed*cf)
    !>   gflow(3)    = its gradient w.r.t. Re, Ncrit top, Ncrit bottom
    !>   kiter       = Newton iterations replayed from the converged state
    subroutine api_adjoint(n, x, y, op, imode, kiter, seed, cp_seed, cf_seed, gx, gy, gop, gflow, res, conv) &
            bind(c, name='xad_adjoint')
        use xad_core_diff, only: xad_eval_b
        use i_xfoil_diff, only: xad_minit, REInf1, ACRit, xad_keepgeom
        use xad_zero_adj, only: xad_zero_adjoints
        use xad_hand_store, only: xad_hand_reset
        integer(c_int), intent(in) :: n, imode, kiter
        real(c_double), intent(in) :: x(n), y(n), op, seed(5), cp_seed(n), cf_seed(n)
        real(c_double), intent(out) :: gx(n), gy(n), gop, gflow(3), res(5)
        integer(c_int), intent(out) :: conv
        real :: xw(n), yw(n), opw, cp(n), cpb(n), cf(n), cfb(n), rew, reb, ncw(2), ncb(2)
        real :: cl, cd, cdf, cm, alfa, clb, cdb, cdfb, cmb, alfab
        logical :: lconv

        xw = x
        yw = y
        opw = op
        rew = REInf1
        ncw = ACRit
        gx = 0.0
        gy = 0.0
        gop = 0.0
        reb = 0.0
        ncb = 0.0
        clb = seed(1)
        cdb = seed(2)
        cdfb = seed(3)
        cmb = seed(4)
        alfab = seed(5)
        cpb = cp_seed
        cfb = cf_seed
        call xad_zero_adjoints()
        call xad_hand_reset()
        ! the reverse sweep restores earlier values of the panel arrays
        geo_ok = .false.
        xad_keepgeom = .false.
        xad_minit = kiter
        call xad_eval_b(n, xw, gx, yw, gy, opw, gop, rew, reb, ncw, ncb, imode, kiter, .true., &
                        cl, clb, cd, cdb, cdf, cdfb, cm, cmb, alfa, alfab, cp, cpb, cf, cfb, lconv)
        xad_minit = 1
        gflow(1) = reb
        gflow(2) = ncb(1)
        gflow(3) = ncb(2)
        res(1) = cl
        res(2) = cd
        res(3) = cdf
        res(4) = cm
        res(5) = alfa
        conv = merge(1, 0, lconv)
    end subroutine api_adjoint


    !> Diagnostic / prototyping: panel and BL data of the last analysis.
    !>   ni(6)  = N, NW, IST, IBLTE(1), IBLTE(2), NBL(2)
    !>   dij(nz, nz) source influence matrix (nz = N + NW <= nzmax)
    !>   xy(nz, 2) node coordinates, qu(nz, 2) QINVU, uedg/uinv/xssi/tau/mass(nbx, 2),
    !>   ipan(nbx, 2), vti(nbx, 2) along the BL of each side, ra(4) = ALFA, CL, CLSPEC, CD
    subroutine api_get_flowdata(nzmax, nbx, ni, dij_o, xy, qu, uedg_o, uinv_o, xssi_o, tau_o, mass_o, ipan_o, &
            vti_o, ra) bind(c, name='xad_get_flowdata')
        use i_xfoil_diff, only: N, NW, IST, IBLte, NBL, DIJ, X, Y, QINvu, UEDg, UINv, XSSi, TAU, MASs, IPAn, VTI, &
                ALFa, CL, CLSpec, CD
        integer(c_int), intent(in) :: nzmax, nbx
        integer(c_int), intent(out) :: ni(6), ipan_o(nbx, 2)
        real(c_double), intent(out) :: dij_o(nzmax, nzmax), xy(nzmax, 2), qu(nzmax, 2), uedg_o(nbx, 2), &
                uinv_o(nbx, 2), xssi_o(nbx, 2), tau_o(nbx, 2), mass_o(nbx, 2), vti_o(nbx, 2), ra(4)
        integer :: nz, k
        nz = N + NW
        ni = (/N, NW, IST, IBLte(1), IBLte(2), NBL(2)/)
        dij_o(1:nz, 1:nz) = DIJ(1:nz, 1:nz)
        xy(1:nz, 1) = X(1:nz)
        xy(1:nz, 2) = Y(1:nz)
        qu(1:nz, :) = QINvu(1:nz, :)
        do k = 1, 2
            uedg_o(1:NBL(k), k) = UEDg(1:NBL(k), k)
            uinv_o(1:NBL(k), k) = UINv(1:NBL(k), k)
            xssi_o(1:NBL(k), k) = XSSi(1:NBL(k), k)
            tau_o(1:NBL(k), k) = TAU(1:NBL(k), k)
            mass_o(1:NBL(k), k) = MASs(1:NBL(k), k)
            ipan_o(1:NBL(k), k) = IPAn(1:NBL(k), k)
            vti_o(1:NBL(k), k) = VTI(1:NBL(k), k)
        enddo
        ra = (/ALFa, CL, CLSpec, CD/)
    end subroutine api_get_flowdata


    !> Size of the BL state buffers of xad_get_state / xad_set_state.
    subroutine api_state_size(nr, ni) bind(c, name='xad_state_size')
        use i_xfoil_diff, only: IVX, ISX, IZX
        integer(c_int), intent(out) :: nr, ni
        nr = 5 * IVX * ISX + IZX + 2
        ni = ISX
    end subroutine api_state_size


    !> Copies the BL state of the last analysis (theta, delta*, ctau / ampl,
    !> Ue, mass defect per station, surface velocities, alpha, CL, transition
    !> intervals), so that a later analysis of a nearby section can restart
    !> from it (xad_set_state, then xad_analyze with warm = 1).
    subroutine api_get_state(rbuf, ibuf) bind(c, name='xad_get_state')
        use i_xfoil_diff, only: IVX, ISX, IZX, THEt, DSTr, CTAu, UEDg, MASs, QVIs, ALFa, CL, ITRan
        real(c_double), intent(out) :: rbuf(5 * IVX * ISX + IZX + 2)
        integer(c_int), intent(out) :: ibuf(ISX)
        integer :: k
        k = IVX * ISX
        rbuf(1:k) = reshape(THEt, (/k/))
        rbuf(k + 1:2 * k) = reshape(DSTr, (/k/))
        rbuf(2 * k + 1:3 * k) = reshape(CTAu, (/k/))
        rbuf(3 * k + 1:4 * k) = reshape(UEDg, (/k/))
        rbuf(4 * k + 1:5 * k) = reshape(MASs, (/k/))
        rbuf(5 * k + 1:5 * k + IZX) = QVIs
        rbuf(5 * k + IZX + 1) = ALFa
        rbuf(5 * k + IZX + 2) = CL
        ibuf = ITRan
    end subroutine api_get_state


    !> Restores a BL state saved by xad_get_state and marks the BL as
    !> initialized: the next xad_analyze with warm = 1 restarts from it.
    subroutine api_set_state(rbuf, ibuf) bind(c, name='xad_set_state')
        use i_xfoil_diff, only: IVX, ISX, IZX, THEt, DSTr, CTAu, UEDg, MASs, QVIs, ALFa, CL, ITRan, LBLini
        real(c_double), intent(in) :: rbuf(5 * IVX * ISX + IZX + 2)
        integer(c_int), intent(in) :: ibuf(ISX)
        integer :: k
        k = IVX * ISX
        THEt = reshape(rbuf(1:k), (/IVX, ISX/))
        DSTr = reshape(rbuf(k + 1:2 * k), (/IVX, ISX/))
        CTAu = reshape(rbuf(2 * k + 1:3 * k), (/IVX, ISX/))
        UEDg = reshape(rbuf(3 * k + 1:4 * k), (/IVX, ISX/))
        MASs = reshape(rbuf(4 * k + 1:5 * k), (/IVX, ISX/))
        QVIs = rbuf(5 * k + 1:5 * k + IZX)
        ALFa = rbuf(5 * k + IZX + 1)
        CL = rbuf(5 * k + IZX + 2)
        ITRan = ibuf
        LBLini = .true.
    end subroutine api_set_state


    !> Library limits: max panel nodes, max BL stations per side.
    subroutine api_limits(nmax, nblmax) bind(c, name='xad_limits')
        use i_xfoil_diff, only: IQX, IVX
        integer(c_int), intent(out) :: nmax, nblmax
        nmax = IQX - 6
        nblmax = IVX
    end subroutine api_limits

end module xad_api
