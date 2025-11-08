! DarkEnergyIDE - Interacting Dark Energy model
! Complete implementation of IDE physics from liaocrane/IDECAMB
! Supports CPL parametrization with dark sector coupling

module DarkEnergyIDE
    use DarkEnergyInterface
    use results
    use constants
    use classes
    use IDEtools
    use CoupledFluidModels
    implicit none
    private

    type, extends(TDarkEnergyEqnOfState) :: TDarkEnergyIDE
        ! IDE-specific parameters
        real(dl) :: beta = 0._dl          ! Coupling strength
        integer :: coupling_form = 1       ! 1=Hrde, 2=Hrc, 3=H0rde, 4=H0rc
        integer :: eos_form = 1            ! 1=CPL, 2=HDE, 3=NADE
        integer :: covariant_form = 1      ! 1=uc (dark matter frame), 2=ude (dark energy frame)
        logical :: w_perturb = .true.      ! Evolve perturbations
        logical :: evolve_vc = .false.     ! Evolve dark matter velocity
        class(CAMBdata), pointer, private :: State => null()
    contains
        procedure :: ReadParams => TDarkEnergyIDE_ReadParams
        procedure, nopass :: PythonClass => TDarkEnergyIDE_PythonClass
        procedure, nopass :: SelfPointer => TDarkEnergyIDE_SelfPointer
        procedure :: Init => TDarkEnergyIDE_Init
        procedure :: BackgroundDensityAndPressure => TDarkEnergyIDE_BackgroundDensityAndPressure
        procedure :: PerturbedStressEnergy => TDarkEnergyIDE_PerturbedStressEnergy
        procedure :: PerturbationEvolve => TDarkEnergyIDE_PerturbationEvolve
        procedure :: PrintFeedback => TDarkEnergyIDE_PrintFeedback
        procedure, private :: IDEout  ! Main interface to IDE physics
    end type TDarkEnergyIDE

    public TDarkEnergyIDE

contains

    subroutine TDarkEnergyIDE_ReadParams(this, Ini)
        use IniObjects
        class(TDarkEnergyIDE) :: this
        class(TIniFile), intent(in) :: Ini

        ! Read base class parameters (w, wa, cs2)
        call this%TDarkEnergyEqnOfState%ReadParams(Ini)

        ! Read IDE-specific parameters
        this%beta = Ini%Read_Double('beta', 0.d0)
        this%coupling_form = Ini%Read_Int('coupling_form', 1)
        this%eos_form = Ini%Read_Int('eos_form', 1)
        this%covariant_form = Ini%Read_Int('covariant_form', 1)

    end subroutine TDarkEnergyIDE_ReadParams

    function TDarkEnergyIDE_PythonClass()
        character(LEN=:), allocatable :: TDarkEnergyIDE_PythonClass
        TDarkEnergyIDE_PythonClass = 'DarkEnergyIDE'
    end function TDarkEnergyIDE_PythonClass

    subroutine TDarkEnergyIDE_SelfPointer(cptr, P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type(TDarkEnergyIDE), pointer :: PType
        class(TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TDarkEnergyIDE_SelfPointer

    ! Setter function for IDE parameters (callable from Python)
    subroutine TDarkEnergyIDE_SetParams(this, beta, coupling_form, eos_form, covariant_form) bind(C)
        use iso_c_binding
        type(c_ptr), value :: this
        real(c_double), value :: beta
        integer(c_int), value :: coupling_form, eos_form, covariant_form
        type(TDarkEnergyIDE), pointer :: ide

        call c_f_pointer(this, ide)
        ide%beta = beta
        ide%coupling_form = coupling_form
        ide%eos_form = eos_form
        ide%covariant_form = covariant_form
    end subroutine TDarkEnergyIDE_SetParams

    subroutine TDarkEnergyIDE_Init(this, State)
        class(TDarkEnergyIDE), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State

        ! Store state reference
        select type(State)
        class is (CAMBdata)
            this%State => State
        end select

        ! Set up CoupledFluidModels configuration
        CFT%WForm = this%eos_form
        CFT%QForm = this%coupling_form
        CFT%CovQForm = this%covariant_form

        ! Set parameters
        CFP%w0 = this%w_lam
        CFP%w1 = this%wa
        CFP%beta = this%beta

        ! Determine if we need to evolve perturbations and dark matter velocity
        call LogicalPertur_CF(this%w_perturb, this%evolve_vc)

        ! Check if this is cosmological constant
        if (this%w_lam == -1._dl .and. this%wa == 0._dl .and. this%beta == 0._dl) then
            this%is_cosmological_constant = .true.
            this%num_perturb_equations = 0
        else
            this%is_cosmological_constant = .false.
            ! Solve background evolution equations
            if (this%w_perturb) then
                call Solve_CF(State)
            end if
            this%num_perturb_equations = 2
        end if

    end subroutine TDarkEnergyIDE_Init

    subroutine TDarkEnergyIDE_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
        class(TDarkEnergyIDE), intent(inout) :: this
        real(dl), intent(in) :: grhov, a
        real(dl), intent(out) :: grhov_t
        real(dl), optional, intent(out) :: w
        real(dl) :: grhoc_t, wde

        if (this%is_cosmological_constant) then
            ! Cosmological constant
            grhov_t = grhov * a * a
            if (present(w)) w = -1_dl
        else
            ! Get densities and EoS from IDE solver
            call this%IDEout(a, grhov_t, grhoc_t, wde=wde)
            if (present(w)) w = wde
        end if

    end subroutine TDarkEnergyIDE_BackgroundDensityAndPressure

    ! Main interface function: get background and perturbation quantities
    subroutine IDEout(this, a, grhov_t, grhoc_t, wde, gQ, gC, gD, ca2)
        class(TDarkEnergyIDE), intent(in) :: this
        real(dl), intent(in) :: a
        real(dl), intent(out) :: grhov_t, grhoc_t
        real(dl), intent(out), optional :: wde, ca2, gQ, gC(3), gD(2)
        real(dl) :: a2, adotoa

        a2 = a*a

        if (.not. this%w_perturb) then
            ! Cosmological constant case
            grhov_t = this%State%grhov * a2
            grhoc_t = this%State%grhoc / a
            if (present(wde)) wde = -1._dl
            if (present(gQ)) gQ = 0._dl
            if (present(gC)) gC = 0._dl
            if (present(gD)) gD = 0._dl
            if (present(ca2)) ca2 = -1._dl
            return
        end if

        ! Get background densities from interpolation table
        call Get_grhodea2_grhoca2(a, grhov_t, grhoc_t)

        ! Compute Hubble parameter
        adotoa = ConformalH(this%State, a, grhov_t, grhoc_t)

        ! Compute equation of state
        if (present(wde)) then
            wde = EoS_CF(a, adotoa, grhov_t, grhoc_t)
        end if

        ! Compute coupling term
        if (present(gQ)) then
            gQ = Coup_CF(this%State, a, adotoa, grhov_t, grhoc_t)
        end if

        ! Compute perturbation coupling terms
        if (present(gC)) then
            call PerturCoupC_CF(this%State, a, adotoa, grhov_t, grhoc_t, gC)
        end if

        if (present(gD)) then
            call PerturCoupD_CF(this%State, a, adotoa, grhov_t, grhoc_t, gD)
        end if

        ! Compute effective sound speed (simplified for constant w)
        if (present(ca2)) then
            ca2 = EoS_CF(a, adotoa, grhov_t, grhoc_t)
        end if

    end subroutine IDEout

    subroutine TDarkEnergyIDE_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        class(TDarkEnergyIDE), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix
        real(dl) :: vT, ca2, grhoc_t, wde, grhov_t_local

        ! Get dark matter and dark energy densities
        call this%IDEout(a, grhov_t_local, grhoc_t)

        if (k > 0.06_dl*adotoa) then
            vT = ay(w_ix)
            dgrhoe = ay(w_ix + 1)
            dgqe = k*vT*grhov_t_local
        else
            ca2 = this%cs2_lam
            dgrhoe = dgrho*grhov_t_local/(gpres_noDE + grhov_t_local)
            dgqe = 3*adotoa*(ca2 - w)*dgrhoe + k*dgq*grhov_t_local/(gpres_noDE + grhov_t_local)
        end if

    end subroutine TDarkEnergyIDE_PerturbedStressEnergy

    subroutine TDarkEnergyIDE_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
        class(TDarkEnergyIDE), intent(in) :: this
        real(dl), intent(inout) :: ayprime(:)
        real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
        integer, intent(in) :: w_ix
        real(dl) :: vT, sigma, dgrhoe, ca2, clxde, vde
        real(dl) :: grhov_t, grhoc_t, wde, gQ, gC(3), gD(2)
        real(dl) :: k2, a2

        a2 = a*a
        k2 = k*k

        ! Get IDE quantities
        call this%IDEout(a, grhov_t, grhoc_t, wde=wde, gQ=gQ, gC=gC, gD=gD, ca2=ca2)

        vT = y(w_ix)
        dgrhoe = y(w_ix + 1)
        clxde = dgrhoe / grhov_t
        vde = vT

        ! Dark energy velocity evolution
        ayprime(w_ix) = -adotoa*vT*(1 - 3*ca2) - k*ca2*dgrhoe/grhov_t - &
            k*(1 + wde)*z + gD(1)*vde/adotoa - gD(2)*vde/adotoa

        ! Dark energy density perturbation evolution
        ayprime(w_ix + 1) = -3*adotoa*(ca2 - wde)*dgrhoe - (1 + wde)*k*grhov_t*vT + &
            gC(1)*clxde - gC(3)*clxde

    end subroutine TDarkEnergyIDE_PerturbationEvolve

    subroutine TDarkEnergyIDE_PrintFeedback(this, FeedbackLevel)
        class(TDarkEnergyIDE) :: this
        integer, intent(in) :: FeedbackLevel

        if (FeedbackLevel > 0) then
            write(*,'("  Dark Energy: IDE model")')
            write(*,'("    w0 = ", F8.4)') this%w_lam
            write(*,'("    wa = ", F8.4)') this%wa
            write(*,'("    beta = ", F8.4)') this%beta
            write(*,'("    coupling_form = ", I2)') this%coupling_form
            write(*,'("    evolve_vc = ", L1)') this%evolve_vc
        end if

    end subroutine TDarkEnergyIDE_PrintFeedback

end module DarkEnergyIDE
