    module DarkEnergyIDE
    use DarkEnergyPPF
    use classes
    use interpolation
    use MassiveNu, only: ThermalNuBack
    implicit none

    private

    integer, parameter :: ide_qform_hrho_de = 1
    integer, parameter :: ide_qform_hrho_c = 2
    integer, parameter :: ide_qform_h0rho_de = 3
    integer, parameter :: ide_qform_h0rho_c = 4

    type, extends(TDarkEnergyPPF) :: TInteractingDarkEnergy
        real(dl) :: w1 = 0._dl
        real(dl) :: beta = 0._dl
        integer :: q_form = ide_qform_hrho_de
        integer :: cov_q_form = 1
        real(dl) :: amin_background = 1e-8_dl
        integer :: n_background = 2000
        logical :: has_background_solution = .false.
        real(dl), private :: grhov0 = 0._dl
        real(dl), private :: grhoc0 = 0._dl
        type(TCubicSpline), private :: de_a4
        type(TCubicSpline), private :: cdm_a4
    contains
    procedure :: ReadParams => TInteractingDarkEnergy_ReadParams
    procedure, nopass :: PythonClass => TInteractingDarkEnergy_PythonClass
    procedure, nopass :: SelfPointer => TInteractingDarkEnergy_SelfPointer
    procedure :: Init => TInteractingDarkEnergy_Init
    procedure :: w_de => TInteractingDarkEnergy_w_de
    procedure :: grho_de => TInteractingDarkEnergy_grho_de
    procedure :: BackgroundCDMDensity => TInteractingDarkEnergy_BackgroundCDMDensity
    procedure, private :: SolveBackground => TInteractingDarkEnergy_SolveBackground
    procedure, private :: BackgroundDerivs => TInteractingDarkEnergy_BackgroundDerivs
    end type TInteractingDarkEnergy

    public TInteractingDarkEnergy
    contains

    subroutine TInteractingDarkEnergy_ReadParams(this, Ini)
    use IniObjects
    class(TInteractingDarkEnergy) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyPPF%ReadParams(Ini)

    this%w1 = this%wa
    if (Ini%HasKey("w1")) this%w1 = Ini%Read_Double("w1")
    if (Ini%HasKey("ide_w1")) this%w1 = Ini%Read_Double("ide_w1")
    this%wa = this%w1

    this%beta = Ini%Read_Double("ide_beta", 0._dl)
    if (Ini%HasKey("beta_cf")) this%beta = Ini%Read_Double("beta_cf")

    this%q_form = Ini%Read_Int("ide_qform", ide_qform_hrho_de)
    if (Ini%HasKey("QForm_CF")) this%q_form = Ini%Read_Int("QForm_CF")

    this%cov_q_form = Ini%Read_Int("ide_covqform", 1)
    if (Ini%HasKey("CovQForm_CF")) this%cov_q_form = Ini%Read_Int("CovQForm_CF")

    call Ini%Read("ide_amin_background", this%amin_background)
    call Ini%Read("ide_n_background", this%n_background)

    end subroutine TInteractingDarkEnergy_ReadParams


    function TInteractingDarkEnergy_PythonClass()
    character(LEN=:), allocatable :: TInteractingDarkEnergy_PythonClass

    TInteractingDarkEnergy_PythonClass = "InteractingDarkEnergy"
    end function TInteractingDarkEnergy_PythonClass


    subroutine TInteractingDarkEnergy_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TInteractingDarkEnergy), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TInteractingDarkEnergy_SelfPointer


    subroutine TInteractingDarkEnergy_Init(this, State)
    use config
    class(TInteractingDarkEnergy), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    if (this%use_tabulated_w) then
        call GlobalError("InteractingDarkEnergy does not support use_tabulated_w", error_unsupported_params)
        return
    end if
    if (this%n_background < 32) this%n_background = 32
    if (this%amin_background <= 0._dl .or. this%amin_background >= 1._dl) then
        call GlobalError("ide_amin_background must satisfy 0 < ide_amin_background < 1", error_unsupported_params)
        return
    end if
    if (this%q_form < ide_qform_hrho_de .or. this%q_form > ide_qform_h0rho_c) then
        call GlobalError("ide_qform must be 1, 2, 3 or 4", error_unsupported_params)
        return
    end if
    if (this%cov_q_form < 1 .or. this%cov_q_form > 2) then
        call GlobalError("ide_covqform must be 1 or 2", error_unsupported_params)
        return
    end if

    this%wa = this%w1
    call this%TDarkEnergyPPF%Init(State)
    this%grhov0 = State%grhov
    this%grhoc0 = State%grhoc

    call this%SolveBackground(State)
    this%has_background_solution = .true.

    if (abs(this%beta) > 1e-14_dl) then
        this%is_cosmological_constant = .false.
        if (.not. this%no_perturbations) this%num_perturb_equations = 1
    end if

    end subroutine TInteractingDarkEnergy_Init


    subroutine TInteractingDarkEnergy_SolveBackground(this, State)
    class(TInteractingDarkEnergy), intent(inout) :: this
    class(TCAMBdata), intent(in) :: State
    integer :: i, n
    real(dl) :: dloga, a_prev, a_now
    real(dl) :: k1_de, k2_de, k3_de, k4_de
    real(dl) :: k1_c, k2_c, k3_c, k4_c
    real(dl), allocatable :: a_desc(:), de_desc(:), c_desc(:)
    real(dl), allocatable :: a_asc(:), de_asc(:), c_asc(:)

    n = this%n_background
    allocate(a_desc(n), de_desc(n), c_desc(n))
    allocate(a_asc(n), de_asc(n), c_asc(n))

    dloga = log(this%amin_background) / real(n - 1, dl)
    a_desc(1) = 1._dl
    de_desc(1) = max(this%grhov0, 1e-40_dl)
    c_desc(1) = max(this%grhoc0, 1e-40_dl)

    do i = 2, n
        a_prev = a_desc(i - 1)
        a_now = exp(real(i - 1, dl) * dloga)
        a_desc(i) = a_now

        call this%BackgroundDerivs(State, a_prev, de_desc(i - 1), c_desc(i - 1), k1_de, k1_c)
        call this%BackgroundDerivs( &
            State, &
            exp(log(a_prev) + dloga / 2._dl), &
            de_desc(i - 1) + dloga * k1_de / 2._dl, &
            c_desc(i - 1) + dloga * k1_c / 2._dl, &
            k2_de, &
            k2_c &
        )
        call this%BackgroundDerivs( &
            State, &
            exp(log(a_prev) + dloga / 2._dl), &
            de_desc(i - 1) + dloga * k2_de / 2._dl, &
            c_desc(i - 1) + dloga * k2_c / 2._dl, &
            k3_de, &
            k3_c &
        )
        call this%BackgroundDerivs( &
            State, &
            a_now, &
            de_desc(i - 1) + dloga * k3_de, &
            c_desc(i - 1) + dloga * k3_c, &
            k4_de, &
            k4_c &
        )

        de_desc(i) = de_desc(i - 1) + dloga * (k1_de + 2._dl * k2_de + 2._dl * k3_de + k4_de) / 6._dl
        c_desc(i) = c_desc(i - 1) + dloga * (k1_c + 2._dl * k2_c + 2._dl * k3_c + k4_c) / 6._dl

        de_desc(i) = max(de_desc(i), 1e-40_dl)
        c_desc(i) = max(c_desc(i), 1e-40_dl)
    end do

    do i = 1, n
        a_asc(i) = a_desc(n - i + 1)
        de_asc(i) = de_desc(n - i + 1)
        c_asc(i) = c_desc(n - i + 1)
    end do

    call this%de_a4%Init(log(a_asc), de_asc)
    call this%cdm_a4%Init(log(a_asc), c_asc)

    deallocate(a_desc, de_desc, c_desc, a_asc, de_asc, c_asc)

    end subroutine TInteractingDarkEnergy_SolveBackground


    subroutine TInteractingDarkEnergy_BackgroundDerivs(this, State, a, de_a4, c_a4, dlog_de_a4, dlog_c_a4)
    class(TInteractingDarkEnergy), intent(in) :: this
    class(TCAMBdata), intent(in) :: State
    real(dl), intent(in) :: a, de_a4, c_a4
    real(dl), intent(out) :: dlog_de_a4, dlog_c_a4
    real(dl) :: a2, grhov_t, grhoc_t, wde, adotoa, gq
    real(dl) :: rho_a4_non_dark, rhonu
    integer :: nu_i

    a2 = a * a
    grhov_t = de_a4 / a2
    grhoc_t = c_a4 / a2

    rho_a4_non_dark = State%grhok * a2 + State%grhob * a + State%grhog + State%grhornomass
    if (State%CP%Num_Nu_massive /= 0) then
        do nu_i = 1, State%CP%Nu_mass_eigenstates
            call ThermalNuBack%rho(a * State%nu_masses(nu_i), rhonu)
            rho_a4_non_dark = rho_a4_non_dark + rhonu * State%grhormass(nu_i)
        end do
    end if
    adotoa = sqrt(max((rho_a4_non_dark + de_a4 + c_a4) / 3._dl, 1e-40_dl))
    wde = this%w_lam + this%wa * (1._dl - a)

    select case (this%q_form)
    case (ide_qform_hrho_de)
        gq = this%beta * adotoa * grhov_t
    case (ide_qform_hrho_c)
        gq = this%beta * adotoa * grhoc_t
    case (ide_qform_h0rho_de)
        gq = this%beta * sqrt(State%grhocrit / 3._dl) * a * grhov_t
    case (ide_qform_h0rho_c)
        gq = this%beta * sqrt(State%grhocrit / 3._dl) * a * grhoc_t
    case default
        gq = 0._dl
    end select

    dlog_de_a4 = (1._dl - 3._dl * wde) * de_a4 + gq * a2 / adotoa
    dlog_c_a4 = c_a4 - gq * a2 / adotoa

    end subroutine TInteractingDarkEnergy_BackgroundDerivs


    function TInteractingDarkEnergy_w_de(this, a)
    class(TInteractingDarkEnergy) :: this
    real(dl), intent(in) :: a
    real(dl) :: TInteractingDarkEnergy_w_de

    TInteractingDarkEnergy_w_de = this%w_lam + this%wa * (1._dl - a)

    end function TInteractingDarkEnergy_w_de


    function TInteractingDarkEnergy_grho_de(this, a) result(grho_de)
    class(TInteractingDarkEnergy) :: this
    real(dl), intent(in) :: a
    real(dl) :: grho_de, loga, de_a4

    if (.not. this%has_background_solution) then
        grho_de = this%TDarkEnergyEqnOfState%grho_de(a)
        return
    end if

    if (a <= 0._dl) then
        grho_de = 0._dl
        return
    end if
    loga = log(a)
    if (loga <= this%de_a4%X(1)) then
        de_a4 = this%de_a4%F(1)
    elseif (loga >= this%de_a4%X(this%de_a4%n)) then
        de_a4 = this%de_a4%F(this%de_a4%n)
    else
        de_a4 = this%de_a4%Value(loga)
    end if
    grho_de = max(de_a4, 0._dl) / max(this%grhov0, 1e-40_dl)

    end function TInteractingDarkEnergy_grho_de


    subroutine TInteractingDarkEnergy_BackgroundCDMDensity(this, grhoc, a, grhoc_t)
    class(TInteractingDarkEnergy), intent(inout) :: this
    real(dl), intent(in) :: grhoc, a
    real(dl), intent(out) :: grhoc_t
    real(dl) :: loga, c_a4

    if (.not. this%has_background_solution .or. a <= 0._dl) then
        call this%TDarkEnergyModel%BackgroundCDMDensity(grhoc, a, grhoc_t)
        return
    end if
    loga = log(a)
    if (loga <= this%cdm_a4%X(1)) then
        c_a4 = this%cdm_a4%F(1)
    elseif (loga >= this%cdm_a4%X(this%cdm_a4%n)) then
        c_a4 = this%cdm_a4%F(this%cdm_a4%n)
    else
        c_a4 = this%cdm_a4%Value(loga)
    end if
    if (a > 1e-12_dl) then
        grhoc_t = grhoc * (c_a4 / max(this%grhoc0, 1e-40_dl)) / (a * a)
    else
        grhoc_t = 0._dl
    end if

    end subroutine TInteractingDarkEnergy_BackgroundCDMDensity

    end module DarkEnergyIDE
