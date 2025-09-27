    module DarkEnergyComposite
    !===============================================================================
    ! EarlyQuintessencePPF: Composite Dark Energy Model
    !
    ! This module implements a composite dark energy model that combines:
    ! 1. Early Dark Energy (EDE): Quintessence field with axion-like potential
    ! 2. Late-time PPF: Parameterized Post-Friedmann w(a) evolution
    !
    ! Physical motivation:
    ! - EDE component provides extra energy injection at high redshifts (z ~ 1000-10000)
    ! - PPF component handles late-time dark energy evolution with w(a) parameterization
    ! - Composite approach allows addressing Hubble tension while maintaining
    !   consistency with late-time observations
    !
    ! Mathematical framework:
    ! - Total dark energy = ρ_EDE(a) + ρ_PPF(a)
    ! - EDE: Scalar field φ with V(φ) = m²f²[1 - cos(φ/f)]ⁿ
    !        At late times: φ → φ_min, so ρ_EDE(a=1) → 0
    ! - PPF: CPL parameterization w(a) = w₀ + wₐ(1-a) with full Ω_DE today
    !
    ! Authors: [Your institution/group]
    ! Date: 2025
    ! Reference: [If publishing, add relevant citations]
    !===============================================================================
    use precision
    use DarkEnergyInterface
    use DarkEnergyPPF
    use Quintessence
    use results
    use classes
    use config
    implicit none

    private

    ! Composite DE: Early Quintessence (EDE) + Late-time PPF
    !
    ! This type combines two dark energy components:
    ! 1. Early Dark Energy: Quintessence field φ with axion-like potential
    !    V(φ) = m²f²[1 - cos(φ/f)]ⁿ + Λ_cc
    ! 2. PPF component: CPL parameterization w(a) = w₀ + wₐ(1-a)
    !
    ! Key features:
    ! - EDE provides early-time energy injection (z ~ 1000-10000)
    ! - PPF handles late-time dark energy evolution
    ! - Smooth transition between regimes
    ! - Numerically stable across full parameter space
    !
    type, extends(TDarkEnergyModel) :: TEarlyQuintessencePPF
        ! =======================================================================
        ! EXPOSED PARAMETERS - Set via Python interface or INI files
        ! =======================================================================
        
        ! Early Dark Energy (Quintessence) parameters
        real(dl) :: n = 3._dl              ! Potential index: V ∝ [1-cos(φ/f)]ⁿ
        real(dl) :: f = 0.05_dl            ! Decay constant f/M_pl (dimensionless)
        real(dl) :: m = 5d-54              ! Mass parameter in Planck units
        real(dl) :: theta_i = 3.1_dl       ! Initial field value φ_i/f
        real(dl) :: frac_lambda0 = 0._dl   ! Cosmological constant fraction (=0 for composite)
        logical :: use_zc = .true.         ! Use zc parameterization instead of theta_i
        real(dl) :: zc = 3000._dl          ! Critical redshift for peak EDE
        real(dl) :: fde_zc = 0._dl         ! EDE fraction f_EDE at z = z_c
        integer :: npoints = 5000          ! Integration points for EDE evolution
        integer :: min_steps_per_osc = 10  ! Minimum steps per oscillation

        ! Late-time PPF parameters  
        real(dl) :: w_lam = -1._dl         ! CPL parameter: w₀ (equation of state today)
        real(dl) :: wa = 0._dl             ! CPL parameter: wₐ = -dw/da|_{a=1}
        real(dl) :: cs2_lam = 1._dl        ! Effective sound speed squared

        ! =======================================================================
        ! INTERNAL STATE - Managed automatically
        ! =======================================================================
        type(TEarlyQuintessence) :: EDE    ! EDE component object
        type(TDarkEnergyPPF)     :: PPF    ! PPF component object
        real(dl) :: Omega_PPF_today = 0._dl ! PPF dark energy density parameter today
        integer :: w_ix_EDE = 0            ! Perturbation equation index for EDE
        integer :: w_ix_PPF = 0            ! Perturbation equation index for PPF
    class(CAMBdata), pointer, private :: State => null()  ! Pointer to CAMB state
    contains
        procedure :: ReadParams => TEarlyQuintessencePPF_ReadParams
        procedure, nopass :: PythonClass => TEarlyQuintessencePPF_PythonClass
        procedure, nopass :: SelfPointer => TEarlyQuintessencePPF_SelfPointer
        procedure :: Init => TEarlyQuintessencePPF_Init
        procedure :: BackgroundDensityAndPressure => TEarlyQuintessencePPF_BackgroundDensityAndPressure
        procedure :: PerturbedStressEnergy => TEarlyQuintessencePPF_PerturbedStressEnergy
        procedure :: PerturbationEvolve => TEarlyQuintessencePPF_PerturbationEvolve
        procedure :: PerturbationInitial => TEarlyQuintessencePPF_PerturbationInitial
        procedure :: Effective_w_wa => TEarlyQuintessencePPF_Effective_w_wa
        procedure :: PrintFeedback => TEarlyQuintessencePPF_PrintFeedback
    end type TEarlyQuintessencePPF

    public TEarlyQuintessencePPF

    contains

    subroutine TEarlyQuintessencePPF_ReadParams(this, Ini)
    use IniObjects
    class(TEarlyQuintessencePPF) :: this
    class(TIniFile), intent(in) :: Ini

    ! Optional INI reading (most users will set via Python)
    call this%TDarkEnergyModel%ReadParams(Ini)

    call Ini%Read('EQPPF_n', this%n)
    call Ini%Read('EQPPF_f', this%f)
    call Ini%Read('EQPPF_m', this%m)
    call Ini%Read('EQPPF_theta_i', this%theta_i)
    call Ini%Read('EQPPF_use_zc', this%use_zc)
    call Ini%Read('EQPPF_zc', this%zc)
    call Ini%Read('EQPPF_fde_zc', this%fde_zc)
    call Ini%Read('EQPPF_npoints', this%npoints)
    call Ini%Read('EQPPF_min_steps_per_osc', this%min_steps_per_osc)

    call Ini%Read('EQPPF_w', this%w_lam)
    call Ini%Read('EQPPF_wa', this%wa)
    call Ini%Read('EQPPF_cs2', this%cs2_lam)

    end subroutine TEarlyQuintessencePPF_ReadParams


    function TEarlyQuintessencePPF_PythonClass()
    character(LEN=:), allocatable :: TEarlyQuintessencePPF_PythonClass
    TEarlyQuintessencePPF_PythonClass = 'EarlyQuintessencePPF'
    end function TEarlyQuintessencePPF_PythonClass


    subroutine TEarlyQuintessencePPF_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TEarlyQuintessencePPF), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TEarlyQuintessencePPF_SelfPointer


    subroutine TEarlyQuintessencePPF_Init(this, State)
    !===========================================================================
    ! Initialize the composite EarlyQuintessencePPF model
    !
    ! Initialization strategy:
    ! 1. EDE component: Handles early-time physics (z >> 100)
    !    - Field φ evolves from initial value, oscillates, then settles
    !    - At late times (z ≈ 0): φ → φ_min, so ρ_EDE → 0 naturally
    ! 2. PPF component: Handles all late-time dark energy (z ≲ 100)
    !    - Gets full Ω_DE today since EDE → 0
    !    - Provides w(a) = w₀ + wₐ(1-a) evolution
    !
    ! This approach is physically motivated: EDE is a transient early phenomenon
    ! that doesn't contribute significantly to today's dark energy budget.
    !===========================================================================
    class(TEarlyQuintessencePPF), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: Omega_EDE_today

    ! CRITICAL: Set composite model's state pointer first
    select type(State)
    class is (CAMBdata)
        this%State => State
        
        ! Set default values to avoid uninitialized access
        this%is_cosmological_constant = .false.
        this%num_perturb_equations = 3  ! EDE(2) + PPF(1)
        
        ! Configure EDE subcomponent from exposed fields
        this%EDE%n = this%n
        this%EDE%f = this%f
        this%EDE%m = this%m
        this%EDE%theta_i = this%theta_i
        ! CRITICAL FIX: Set frac_lambda0=0 to avoid double counting with PPF
        this%EDE%frac_lambda0 = 0._dl   ! EDE should have NO lambda contribution
        this%EDE%use_zc = this%use_zc
        this%EDE%zc = this%zc
        this%EDE%fde_zc = this%fde_zc
        this%EDE%npoints = this%npoints
        this%EDE%min_steps_per_osc = this%min_steps_per_osc

        ! Initialize EDE component
        call this%EDE%Init(State)

        ! Configure PPF subcomponent
        this%PPF%w_lam = this%w_lam
        this%PPF%wa = this%wa
        this%PPF%cs2_lam = this%cs2_lam
        this%PPF%use_tabulated_w = .false.
        
        ! Initialize PPF component  
        call this%PPF%Init(State)
        
        ! Dark energy density allocation:
        ! EDE naturally decays to zero at late times (φ settles to minimum)
        ! PPF handles all dark energy density today
        Omega_EDE_today = 0._dl                    ! EDE contribution today ≈ 0
        this%Omega_PPF_today = State%Omega_de      ! PPF gets all dark energy

    end select

    ! Decide number of perturbation equations: EDE (2), PPF (1 if not cosmological constant)
    this%w_ix_EDE = 0
    this%w_ix_PPF = 0

    if (this%EDE%is_cosmological_constant) then
        ! EDE as CC is unlikely given frac_lambda0=0, but handle anyway
        this%EDE%num_perturb_equations = 0
    else
        this%EDE%num_perturb_equations = 2
    end if

    if (this%Omega_PPF_today <= 0._dl) then
        this%PPF%num_perturb_equations = 0
    else
        if (abs(this%PPF%w_lam + 1._dl) < 1.e-6_dl .and. this%PPF%wa==0._dl) then
            this%PPF%num_perturb_equations = 0  ! behaves like CC
        else
            this%PPF%num_perturb_equations = 1
        end if
    end if

    this%num_perturb_equations = this%EDE%num_perturb_equations + this%PPF%num_perturb_equations

    ! Assign indices (relative; the caller provides base w_ix)
    if (this%EDE%num_perturb_equations > 0) then
        this%w_ix_EDE = 0
        this%w_ix_PPF = this%EDE%num_perturb_equations
    else
        this%w_ix_EDE = -1
        this%w_ix_PPF = 0
    end if

    this%is_cosmological_constant = .false.
    end subroutine TEarlyQuintessencePPF_Init


    subroutine TEarlyQuintessencePPF_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    class(TEarlyQuintessencePPF), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w
    real(dl) :: grhov_t_EDE, wE, grhov_t_PPF, wP, grho_de_PPF

    ! Safety check
    if (.not. associated(this%State)) then
        grhov_t = 0._dl
        if (present(w)) w = -1._dl
        return
    end if

    ! Early Quintessence contribution
    call this%EDE%BackgroundDensityAndPressure(grhov, a, grhov_t_EDE, wE)

    ! Late-time PPF background: use PPF's built-in grho_de function
    if (this%Omega_PPF_today > 0._dl) then
        ! Use PPF's internal grho_de calculation (handles numerical stability)
        grho_de_PPF = this%PPF%grho_de(a)
        
        ! Convert to grhov_t using PPF's standard formula
        select type(State => this%State)
        class is (CAMBdata)
            if (a > 1.e-10_dl) then
                grhov_t_PPF = this%Omega_PPF_today * State%grhocrit * grho_de_PPF / (a*a)
            else
                grhov_t_PPF = 0._dl
            end if
        end select
        wP = this%PPF%w_de(a)
    else
        grhov_t_PPF = 0._dl
        wP = -1._dl
    end if

    grhov_t = grhov_t_EDE + grhov_t_PPF
    if (present(w)) then
        if (grhov_t > 0) then
            w = (wE*grhov_t_EDE + wP*grhov_t_PPF) / grhov_t
        else
            w = -1._dl
        end if
    end if
    end subroutine TEarlyQuintessencePPF_BackgroundDensityAndPressure


    subroutine TEarlyQuintessencePPF_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TEarlyQuintessencePPF), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) :: dgrhoe_E, dgqe_E, dgrhoe_P, dgqe_P
    real(dl) :: grhov_t_EDE, wE, grhov_t_PPF, wP, fCPL

    dgrhoe = 0._dl
    dgqe   = 0._dl

    ! EDE perturbations: use EDE's actual background density
    if (this%EDE%num_perturb_equations > 0) then
        ! Get EDE's actual background contribution
        select type(State => this%State)
        class is (CAMBdata)
            call this%EDE%BackgroundDensityAndPressure(State%grhov, a, grhov_t_EDE, wE)
        end select
        call this%EDE%PerturbedStressEnergy(dgrhoe_E, dgqe_E, a, dgq, dgrho, grho, grhov_t_EDE, wE, &
            gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix + this%w_ix_EDE)
        dgrhoe = dgrhoe + dgrhoe_E
        dgqe   = dgqe   + dgqe_E
    end if

    ! PPF perturbations: use PPF's actual background density
    if (this%PPF%num_perturb_equations > 0 .and. this%Omega_PPF_today > 0._dl) then
        ! Calculate PPF's actual background contribution (same as BackgroundDensityAndPressure)
        if (a > 1.e-10_dl) then
            grhov_t_PPF = this%Omega_PPF_today * this%State%grhocrit * this%PPF%grho_de(a) / (a*a)
        else
            grhov_t_PPF = 0._dl
        end if
        wP = this%PPF%w_de(a)
        call this%PPF%PerturbedStressEnergy(dgrhoe_P, dgqe_P, a, dgq, dgrho, grho, grhov_t_PPF, wP, &
            gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix + this%w_ix_PPF)
        dgrhoe = dgrhoe + dgrhoe_P
        dgqe   = dgqe   + dgqe_P
    end if
    end subroutine TEarlyQuintessencePPF_PerturbedStressEnergy


    subroutine TEarlyQuintessencePPF_PerturbationEvolve(this, ayprime, w, w_ix, a, adotoa, k, z, y)
    class(TEarlyQuintessencePPF), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: w
    integer, intent(in) :: w_ix
    real(dl), intent(in) :: a, adotoa, k, z
    real(dl), intent(in) :: y(:)

    if (this%EDE%num_perturb_equations > 0) then
        ! use component w for EDE
        call this%EDE%PerturbationEvolve(ayprime, this%EDE%w_de(a), w_ix + this%w_ix_EDE, a, adotoa, k, z, y)
    end if
    if (this%PPF%num_perturb_equations > 0) then
        ! use component w for PPF
        call this%PPF%PerturbationEvolve(ayprime, this%PPF%w_de(a), w_ix + this%w_ix_PPF, a, adotoa, k, z, y)
    end if
    end subroutine TEarlyQuintessencePPF_PerturbationEvolve


    subroutine TEarlyQuintessencePPF_PerturbationInitial(this, y, a, tau, k)
    class(TEarlyQuintessencePPF), intent(in) :: this
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: a, tau, k

    if (this%EDE%num_perturb_equations > 0) then
        call this%EDE%PerturbationInitial(y, a, tau, k)
    end if
    if (this%PPF%num_perturb_equations > 0) then
        call this%PPF%PerturbationInitial(y, a, tau, k)
    end if
    end subroutine TEarlyQuintessencePPF_PerturbationInitial


    subroutine TEarlyQuintessencePPF_Effective_w_wa(this, w, wa)
    class(TEarlyQuintessencePPF), intent(inout) :: this
    real(dl), intent(out) :: w, wa
    real(dl) :: wppf, wappf

    ! Default to late-time CPL effective parameters (EDE small today)
    call this%PPF%Effective_w_wa(wppf, wappf)
    w = wppf
    wa = wappf
    end subroutine TEarlyQuintessencePPF_Effective_w_wa


    subroutine TEarlyQuintessencePPF_PrintFeedback(this, FeedbackLevel)
    class(TEarlyQuintessencePPF) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("EarlyQuintessencePPF: (w₀, wₐ) = (", f8.5, ", ", f8.5, ")  Ω_PPF,0=", f9.6)') &
            this%w_lam, this%wa, this%Omega_PPF_today
        write(*,'("  EDE parameters: n=", f6.2, ", f=", f8.5, ", use_zc=", L1)') &
            this%n, this%f, this%use_zc
        if (this%use_zc) then
            write(*,'("  zc=", f8.0, ", fde_zc=", f8.5)') this%zc, this%fde_zc
        else
            write(*,'("  theta_i=", f8.3, ", m=", es12.3)') this%theta_i, this%m
        end if
    end if
    end subroutine TEarlyQuintessencePPF_PrintFeedback

    end module DarkEnergyComposite
