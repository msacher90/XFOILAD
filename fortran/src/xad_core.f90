!***********************************************************************
!   XFOILAD - adjoint-enabled XFOIL core
!
!   Driver routines around the XFOIL modules (D. de Vries f90 port, GPL).
!   XAD_EVAL is the head routine differentiated by Tapenade in reverse
!   mode: (x_in, y_in, op_in, re_in, nc_in(2)) -> (cl, cd, cdf, cm, alfa, cp(:), cf(:)).
!
!   Gradient strategy ("reverse accumulation from the converged state"):
!     1) XAD_EVAL(warm=.false.) converges the viscous Newton problem from
!        scratch (not differentiated).
!     2) XAD_EVAL_B(warm=.true.) re-runs the whole geometry/panel chain and
!        only K Newton iterations restarted from the converged BL state.
!        Since the Newton map is contractive at the solution, the reverse
!        sweep of these K iterations gives the adjoint gradient of the
!        converged outputs (error ~ |I - J_xfoil^-1 J|^K).
!***********************************************************************
module xad_core
    implicit none
contains


subroutine xad_eval(n_in, x_in, y_in, op_in, re_in, nc_in, imode, niter, warm, &
                    cl_out, cd_out, cdf_out, cm_out, alfa_out, cp_out, cf_out, conv)
    use i_xfoil
    use m_spline, only: scalc, segspl
    use m_xgdes, only: abcopy
    use m_xoper_AD, only: specal, speccl, viscal
    implicit none
    integer, intent(in) :: n_in, imode, niter
    logical, intent(in) :: warm
    real, intent(in) :: x_in(n_in), y_in(n_in), op_in, re_in, nc_in(2)
    real, intent(out) :: cl_out, cd_out, cdf_out, cm_out, alfa_out
    real, intent(out) :: cp_out(n_in), cf_out(n_in)
    logical, intent(out) :: conv
    integer :: i, is, ibl
    !
    !---- Reynolds number and critical amplification factors (top, bottom):
    !     inputs of the head routine, so that the adjoint gives their gradients
    REInf1 = re_in
    ACRit(1) = nc_in(1)
    ACRit(2) = nc_in(2)
    !
    !---- buffer airfoil = panel nodes (no repaneling: nodes come from caller,
    !     counterclockwise from upper TE around LE to lower TE)
    if (xad_keepgeom) then
        !----- same nodes as the previous analysis: geometry, factored AIJ, GAMU
        !      and airfoil block of DIJ kept, everything else reset as in ABCOPY
        LQInu = .false.
        LWAke = .false.
        LWDij = .false.
        LIPan = .false.
        LVConv = .false.
        LSCini = .false.
    else
        do i = 1, n_in
            XB(i) = x_in(i)
            YB(i) = y_in(i)
        enddo
        NB = n_in
        call scalc(XB, YB, SB, NB)
        call segspl(XB, XBP, SB, NB)
        call segspl(YB, YBP, SB, NB)
        !
        !---- current airfoil from buffer, resets all geometry-dependent flags
        call abcopy(.false.)
    endif
    !
    if (.not.warm) then
        LBLini = .false.
    endif
    !---- forget the transition-cycle under-relaxation of the previous call
    do i = 1, 2
        xad_trprv(i) = 0
        xad_trnfl(i) = 0
        xad_mfprv(i) = 0
        xad_mfcnt(i) = 0
        xad_mfcur(i) = 0
    enddo
    xad_rlxcap = 1.0
    xad_trstab = 0
    LIPan = .false.
    LVConv = .false.
    QINf = 1.0
    !
    if (imode==1) then
        !----- specified alpha (deg)
        ADEg = op_in
        ALFa = op_in * DTOr
        LALfa = .true.
        call specal
    else
        !----- specified CL (alpha is the Newton unknown)
        CLSpec = op_in
        !------ inviscid alpha for CL: XFOIL builds the wake at this alpha and
        !       keeps it frozen during the viscous Newton iterations, so the
        !       warm restart must do the same to differentiate the same function
        call speccl
        LALfa = .false.
    endif
    !
    if (LVIsc) then
        conv = viscal(niter)
        conv = LVConv .and. conv
    else
        conv = .true.
    endif
    !
    cl_out = CL
    cd_out = CD
    cdf_out = CDF
    cm_out = CM
    alfa_out = ALFa / DTOr
    if (LVIsc) then
        do i = 1, n_in
            cp_out(i) = CPV(i)
        enddo
    else
        do i = 1, n_in
            cp_out(i) = CPI(i)
        enddo
    endif
    !
    !---- skin friction coefficient at the panel nodes, tau_w / (q_inf), signed
    !     along the BL marching direction (< 0: separated flow).  Every node
    !     is a BL station of one side: IPAN(2, 1) = IST, IPAN(2, 2) = IST + 1
    do i = 1, n_in
        cf_out(i) = 0.0
    enddo
    if (LVIsc) then
        do is = 1, 2
            do ibl = 2, IBLte(is)
                cf_out(IPAn(ibl, is)) = TAU(ibl, is) * 2.0 / QINf**2
            enddo
        enddo
    endif
end subroutine xad_eval



subroutine xad_init()
    use s_xfoil, only: mrcl, comset
    use m_xbl, only: blpini, preptrs
    use i_xfoil
    implicit none
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Local variables
    !
    real :: ann
    integer :: i, ip, k, nn
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
    !---- prepare BL pointers
    call preptrs
    !
    !---------------------------------------------------
    !     Variable initialization/default routine.
    !     See file XFOIL.INC for variable description.
    !---------------------------------------------------
    !
    PI = 4.0 * atan(1.0)
    HOPi = 0.50 / PI
    QOPi = 0.25 / PI
    DTOr = PI / 180.0
    !
    !---- default Cp/Cv (air)
    GAMma = 1.4
    GAMm1 = GAMma - 1.0
    !
    !---- set unity freestream speed
    QINf = 1.0
    !
    !---- initialize freestream Mach number to zero
    MATyp = 1
    MINf1 = 0.
    !
    ALFa = 0.0
    COSa = 1.0
    SINa = 0.0
    !
    do i = 1, IQX
        GAMu(i, 1) = 0.
        GAMu(i, 2) = 0.
        GAM(i) = 0.
        GAM_a(i) = 0.
    enddo
    PSIo = 0.
    !
    CL = 0.
    CM = 0.
    CD = 0.
    !
    SIGte = 0.0
    GAMte = 0.0
    SIGte_a = 0.
    GAMte_a = 0.
    !
    do i = 1, IZX
        SIG(i) = 0.
    enddo
    !
    NQSp = 0
    do k = 1, IPX
        ALQsp(k) = 0.
        CLQsp(k) = 0.
        CMQsp(k) = 0.
        do i = 1, IBX
            QSPec(i, k) = 0.
        enddo
    enddo
    !
    AWAke = 0.0
    AVIsc = 0.0
    !
    KIMage = 1
    YIMage = -10.0
    LIMage = .false.
    !
    !---- output coordinate file delimiters
    !-    KDELIM =  0  blanks
    !-           =  1  commas
    !-           =  2  tabs
    KDElim = 0
    !
    LGAmu = .false.
    LQInu = .false.
    LVIsc = .false.
    LWAke = .false.
    LPAcc = .false.
    LBLini = .false.
    LIPan = .false.
    LQAij = .false.
    LADij = .false.
    LWDij = .false.
    LCPxx = .false.
    LQVdes = .false.
    LQSpec = .false.
    LQRefl = .false.
    LVConv = .false.
    LCPref = .false.
    LFOref = .false.
    LPFile = .false.
    LPFilx = .false.
    LPPsho = .false.
    LBFlap = .false.
    LFLap = .false.
    LEIw = .false.
    LSCini = .false.
    LPLot = .false.
    LCLip = .false.
    LVLab = .true.
    LCMinp = .false.
    LHMomp = .false.
    LFReqp = .true.
    !
    LCUrs = .true.
    LLAnd = .true.
    LGSame = .false.
    !
    LGParm = .true.
    LPLcam = .false.
    !
    !---- input airfoil will not be normalized
    LNOrm = .false.
    !
    !---- airfoil will not be forced symmetric
    LQSym = .false.
    LGSym = .false.
    !
    !---- endpoint slopes will be matched
    LQSlop = .true.
    LGSlop = .true.
    LCSlop = .true.
    !
    !---- grids on Qspec(s) and buffer airfoil geometry plots will be plotted
    LQGrid = .true.
    LGGrid = .true.
    LGTick = .true.
    !
    !---- no grid on Cp plots
    LCPgrd = .false.
    !
    !---- overlay inviscid Cp on viscous Cp(x) plot
    LCPinv = .true.

    !---- grid and no symbols are to be used on BL variable plots
    LBLgrd = .true.
    LBLsym = .false.
    !
    !---- buffer and current airfoil flap hinge coordinates
    XBF = 0.0
    YBF = 0.0
    XOF = 0.0
    YOF = 0.0
    !
    NCPref = 0
    !                                               n
    !---- circle plane array size (257, or largest 2  + 1 that will fit array size)
    ann = log(float((2 * IQX) - 1)) / log(2.0)
    nn = int(ann + 0.00001)
    NC1 = 2**nn + 1
    NC1 = min(NC1, 257)
    !
    !---- default paneling parameters
    NPAn = 160
    CVPar = 1.0
    CTErat = 0.15
    CTRrat = 0.2
    !
    !---- default paneling refinement zone x/c endpoints
    XSRef1 = 1.0
    XSRef2 = 1.0
    XPRef1 = 1.0
    XPRef2 = 1.0
    !
    !---- no polars present to begin with
    NPOl = 0
    IPAct = 0
    do ip = 1, NPX
        PFName(ip) = ' '
        PFNamx(ip) = ' '
    enddo
    !
    !---- no reference polars
    NPOlref = 0
    !
    !---- plot aspect ratio, character size
    PLOtar = 0.55
    CH = 0.015
    !
    !---- airfoil node tick-mark size (as fraction of arc length)
    GTIck = 0.0005
    !
    !---- Cp limits in  Cp vs x  plot
    CPMax = 1.0
    CPMin = -2.0
    CPDel = -0.5
    PFAc = PLOtar / (CPMax - CPMin)
    !
    !---- Ue limits in  Ue vs x  plot
    UEMax = 1.8
    UEMin = -1.0
    UEDel = 0.2
    UFAc = PLOtar / (UEMax - UEMin)
    !
    !---- DCp limits in CAMB loading plot
    YPMin = -0.4
    YPMax = 0.4
    !
    !---- scaling factor for Cp vector plot
    VFAc = 0.25
    !
    !---- offsets and scale factor for airfoil in  Cp vs x  plot
    XOFair = 0.09
    YOFair = -.01
    FACair = 0.70
    !
    !---- u/Qinf scale factor for profile plotting
    UPRwt = 0.02
    !
    !---- polar variables to be written to polar save file
    IPOl(1) = IAL
    IPOl(2) = ICL
    IPOl(3) = ICD
    IPOl(4) = ICP
    IPOl(5) = ICM
    NIPol = 5
    NIPol0 = 5
    !
    JPOl(1) = JTN
    JPOl(2) = JTI
    NJPol = 2
    !
    !---- default Cm reference location
    XCMref = 0.25
    YCMref = 0.
    !
    !---- default viscous parameters
    RETyp = 1
    REInf1 = 0.
    ACRit(1) = 9.0
    ACRit(2) = 9.0
    XSTrip(1) = 1.0
    XSTrip(2) = 1.0
    XOCtr(1) = 1.0
    XOCtr(2) = 1.0
    YOCtr(1) = 0.
    YOCtr(2) = 0.
    WAKlen = 1.0
    !
    IDAmp = 0
    !
    !---- select default BL plotting coordinate (can be changed in VPLO)
    IXBlp = 1       !  x
    !cc   IXBLP = 2   !  s
    !
    !---- set BL calibration parameters
    call blpini
    !
    !---- Newton iteration limit
    ITMax = 20
    !
    !---- max number of unconverged sequence points for early exit
    NSEqex = 4
    !
    !---- drop tolerance for BL system solver
    VACcel = 0.01
    !
    !---- inverse-mapping auto-filter level
    FFIlt = 0.0
    !
    !---- default overlay airfoil filename
    ONAme = ' '
    !
    !---- default filename prefix
    PREfix = ' '
    !
    !
    NNAme = 32
    NAMe = '                                '
    !CC             12345678901234567890123456789012
    !
    !---- MSES domain parameters (not used in XFOIL)
    ISPars = ' -2.0  3.0  -2.5  3.5'
    !
    !---- set MINF, REINF, based on current CL-dependence
    call mrcl(1.0, MINf_cl, REInf_cl)
    !
    !---- set various compressibility parameters from MINF
    call comset
    !
    show_output = .false.
    !
end subroutine xad_init


end module xad_core
