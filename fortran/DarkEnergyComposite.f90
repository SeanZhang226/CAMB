    module DarkEnergyComposite
    use precision
    use DarkEnergyInterface
    use DarkEnergyPPF
    use Quintessence
    use results, only: CAMBdata
    use classes
    use config
    implicit none

    private
    real(dl), parameter :: a_eval_min = 1.e-10_dl
    real(dl), parameter :: w_cosmo_constant_tol = 1.e-6_dl

    type, extends(TDarkEnergyModel) :: TEarlyQuintessencePPF
        ! Early dark energy (quintessence) parameters
        real(dl) :: n = 3._dl
        real(dl) :: f = 0.05_dl
        real(dl) :: m = 5e-54_dl
        real(dl) :: theta_i = 3.1_dl
        real(dl) :: frac_lambda0 = 0._dl
        logical :: use_zc = .true.
        real(dl) :: zc = 3000._dl
        real(dl) :: fde_zc = 0._dl
        integer :: npoints = 5000
        integer :: min_steps_per_osc = 10

        ! Late-time PPF parameters
        real(dl) :: w_lam = -1._dl
        real(dl) :: wa = 0._dl
        real(dl) :: cs2_lam = 1._dl

        ! Internal state
        type(TEarlyQuintessence) :: EDE
        type(TDarkEnergyPPF) :: PPF
        real(dl) :: Omega_PPF_today = 0._dl
        integer :: w_ix_EDE = 0
        integer :: w_ix_PPF = 0
        type(CAMBdata), pointer, private :: State => null()
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

    call this%TDarkEnergyModel%ReadParams(Ini)

    call Ini%Read("EQPPF_n", this%n)
    call Ini%Read("EQPPF_f", this%f)
    call Ini%Read("EQPPF_m", this%m)
    call Ini%Read("EQPPF_theta_i", this%theta_i)
    call Ini%Read("EQPPF_use_zc", this%use_zc)
    call Ini%Read("EQPPF_zc", this%zc)
    call Ini%Read("EQPPF_fde_zc", this%fde_zc)
    call Ini%Read("EQPPF_npoints", this%npoints)
    call Ini%Read("EQPPF_min_steps_per_osc", this%min_steps_per_osc)
    call Ini%Read("EQPPF_w", this%w_lam)
    call Ini%Read("EQPPF_wa", this%wa)
    call Ini%Read("EQPPF_cs2", this%cs2_lam)

    end subroutine TEarlyQuintessencePPF_ReadParams


    function TEarlyQuintessencePPF_PythonClass()
    character(LEN=:), allocatable :: TEarlyQuintessencePPF_PythonClass

    TEarlyQuintessencePPF_PythonClass = "EarlyQuintessencePPF"
    end function TEarlyQuintessencePPF_PythonClass


    subroutine TEarlyQuintessencePPF_SelfPointer(cptr, P)
    use iso_c_binding
    type(c_ptr) :: cptr
    type(TEarlyQuintessencePPF), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TEarlyQuintessencePPF_SelfPointer


    subroutine TEarlyQuintessencePPF_Init(this, State)
    class(TEarlyQuintessencePPF), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%State => null()
    select type(State)
    type is (CAMBdata)
        this%State => State
    class default
        call GlobalError("EarlyQuintessencePPF requires CAMBdata state", error_unsupported_params)
        return
    end select

    this%EDE%n = this%n
    this%EDE%f = this%f
    this%EDE%m = this%m
    this%EDE%theta_i = this%theta_i
    this%EDE%frac_lambda0 = this%frac_lambda0
    this%EDE%use_zc = this%use_zc
    this%EDE%zc = this%zc
    this%EDE%fde_zc = this%fde_zc
    this%EDE%npoints = this%npoints
    this%EDE%min_steps_per_osc = this%min_steps_per_osc
    call this%EDE%Init(State)

    this%PPF%w_lam = this%w_lam
    this%PPF%wa = this%wa
    this%PPF%cs2_lam = this%cs2_lam
    this%PPF%use_tabulated_w = .false.
    call this%PPF%Init(State)

    ! EDE decays away at late times; PPF captures late-time DE budget.
    this%Omega_PPF_today = max(this%State%Omega_de - this%frac_lambda0, 0._dl)

    if (this%EDE%is_cosmological_constant) then
        this%EDE%num_perturb_equations = 0
    else
        this%EDE%num_perturb_equations = 2
    end if

    if (this%Omega_PPF_today <= 0._dl) then
        this%PPF%num_perturb_equations = 0
    else
        if (abs(this%PPF%w_lam + 1._dl) < w_cosmo_constant_tol .and. this%PPF%wa == 0._dl) then
            this%PPF%num_perturb_equations = 0
        else
            this%PPF%num_perturb_equations = 1
        end if
    end if

    this%num_perturb_equations = this%EDE%num_perturb_equations + this%PPF%num_perturb_equations

    this%w_ix_EDE = 0
    this%w_ix_PPF = this%EDE%num_perturb_equations
    this%is_cosmological_constant = .false.

    end subroutine TEarlyQuintessencePPF_Init


    subroutine TEarlyQuintessencePPF_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    class(TEarlyQuintessencePPF), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w
    real(dl) :: grhov_t_EDE, wE, grhov_t_PPF, wP

    if (.not. associated(this%State)) then
        grhov_t = 0._dl
        if (present(w)) w = -1._dl
        return
    end if

    call this%EDE%BackgroundDensityAndPressure(grhov, a, grhov_t_EDE, wE)

    if (this%Omega_PPF_today > 0._dl .and. a > a_eval_min) then
        grhov_t_PPF = this%Omega_PPF_today * this%State%grhocrit * this%PPF%grho_de(a) / (a * a)
        wP = this%PPF%w_de(a)
    else
        grhov_t_PPF = 0._dl
        wP = -1._dl
    end if

    grhov_t = grhov_t_EDE + grhov_t_PPF
    if (present(w)) then
        if (grhov_t > 0._dl) then
            w = (wE * grhov_t_EDE + wP * grhov_t_PPF) / grhov_t
        else
            w = -1._dl
        end if
    end if

    end subroutine TEarlyQuintessencePPF_BackgroundDensityAndPressure


    subroutine TEarlyQuintessencePPF_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TEarlyQuintessencePPF), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) :: dgrhoe_E, dgqe_E, dgrhoe_P, dgqe_P
    real(dl) :: grhov_t_EDE, wE, grhov_t_PPF, wP

    dgrhoe = 0._dl
    dgqe = 0._dl

    if (this%EDE%num_perturb_equations > 0) then
        call this%EDE%BackgroundDensityAndPressure(this%State%grhov, a, grhov_t_EDE, wE)
        call this%EDE%PerturbedStressEnergy(dgrhoe_E, dgqe_E, a, dgq, dgrho, grho, grhov_t_EDE, wE, &
            gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix + this%w_ix_EDE)
        dgrhoe = dgrhoe + dgrhoe_E
        dgqe = dgqe + dgqe_E
    end if

    if (this%PPF%num_perturb_equations > 0 .and. this%Omega_PPF_today > 0._dl) then
        if (a > a_eval_min) then
            grhov_t_PPF = this%Omega_PPF_today * this%State%grhocrit * this%PPF%grho_de(a) / (a * a)
        else
            grhov_t_PPF = 0._dl
        end if
        wP = this%PPF%w_de(a)
        call this%PPF%PerturbedStressEnergy(dgrhoe_P, dgqe_P, a, dgq, dgrho, grho, grhov_t_PPF, wP, &
            gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix + this%w_ix_PPF)
        dgrhoe = dgrhoe + dgrhoe_P
        dgqe = dgqe + dgqe_P
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
        call this%EDE%PerturbationEvolve(ayprime, this%EDE%w_de(a), w_ix + this%w_ix_EDE, a, adotoa, k, z, y)
    end if
    if (this%PPF%num_perturb_equations > 0) then
        call this%PPF%PerturbationEvolve(ayprime, this%PPF%w_de(a), w_ix + this%w_ix_PPF, a, adotoa, k, z, y)
    end if

    end subroutine TEarlyQuintessencePPF_PerturbationEvolve


    subroutine TEarlyQuintessencePPF_PerturbationInitial(this, y, a, tau, k)
    class(TEarlyQuintessencePPF), intent(in) :: this
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer :: i0

    if (this%EDE%num_perturb_equations > 0) then
        call this%EDE%PerturbationInitial(y(1:this%EDE%num_perturb_equations), a, tau, k)
    end if
    if (this%PPF%num_perturb_equations > 0) then
        i0 = this%EDE%num_perturb_equations + 1
        call this%PPF%PerturbationInitial(y(i0:i0 + this%PPF%num_perturb_equations - 1), a, tau, k)
    end if

    end subroutine TEarlyQuintessencePPF_PerturbationInitial


    subroutine TEarlyQuintessencePPF_Effective_w_wa(this, w, wa)
    class(TEarlyQuintessencePPF), intent(inout) :: this
    real(dl), intent(out) :: w, wa

    call this%PPF%Effective_w_wa(w, wa)

    end subroutine TEarlyQuintessencePPF_Effective_w_wa


    subroutine TEarlyQuintessencePPF_PrintFeedback(this, FeedbackLevel)
    class(TEarlyQuintessencePPF) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*, '("EarlyQuintessencePPF: (w0, wa)=(", f8.5, ", ", f8.5, ") Omega_PPF,0=", f9.6)') &
            this%w_lam, this%wa, this%Omega_PPF_today
        write(*, '("  EDE parameters: n=", f6.2, ", f=", f8.5, ", use_zc=", L1)') this%n, this%f, this%use_zc
        if (this%use_zc) then
            write(*, '("  zc=", f8.0, ", fde_zc=", f8.5)') this%zc, this%fde_zc
        else
            write(*, '("  theta_i=", f8.3, ", m=", es12.3)') this%theta_i, this%m
        end if
    end if

    end subroutine TEarlyQuintessencePPF_PrintFeedback

    end module DarkEnergyComposite
