! CoupledFluidModels module for calculating background evolution of
! coupled fluid (Interacting Dark Energy) models
! Ported from liaocrane/IDECAMB

module CoupledFluidModels
    use precision
    use IDEtools
    use classes
    use results
    implicit none

    ! Model types
    integer, parameter :: CPL = 1, HDE = 2, NADE = 3  ! Dark energy EoS forms
    integer, parameter :: Hrde = 1, Hrc = 2, H0rde = 3, H0rc = 4  ! Coupling forms
    integer, parameter :: uc = 1, ude = 2  ! Covariant forms

    ! Model configuration type
    type CoupFluidTypes
        integer :: WForm = 1      ! EoS form (CPL=1, HDE=2, NADE=3)
        integer :: QForm = 1      ! Coupling form (Hrde=1, Hrc=2, H0rde=3, H0rc=4)
        integer :: CovQForm = 1   ! Covariant form (uc=1, ude=2)
    end type CoupFluidTypes
    type(CoupFluidTypes) :: CFT

    ! Model parameters
    type CoupFluidParams
        real(dl) :: w0 = -1_dl    ! CPL model parameter w_0
        real(dl) :: w1 = 0._dl    ! CPL model parameter w_a (using w1 to match IDECAMB)
        real(dl) :: c = 0.7_dl    ! HDE model parameter
        real(dl) :: beta = 0._dl  ! Dimensionless coupling constant
    end type CoupFluidParams
    type(CoupFluidParams) :: CFP

    ! Iteration parameters (for NADE model)
    type IterationParams
        real(dl) :: n = 3.0_dl    ! NADE model parameter
        real(dl) :: diff = 0._dl
    end type IterationParams
    type(IterationParams), save :: ItP

    ! Background evolution arrays
    integer, parameter :: nvar = 2   ! Number of variables: rho_de, rho_c
    integer, parameter :: na = 2000  ! Number of time steps
    integer, parameter :: i_de4 = 1, i_c4 = 2  ! Array indices
    real(dl), parameter :: amin = 1.d-12  ! Minimum scale factor

    real(dl) :: ai(na)        ! Scale factors a_i
    real(dl) :: ya(na, nvar)  ! Function values [grhov_t*a^2, grhoc_t*a^2]
    real(dl) :: dya(na, nvar) ! Derivatives d/da of ya
    real(dl) :: a_trans       ! Transition scale factor (when rho_de goes negative)

    ! State reference for accessing cosmological parameters
    class(CAMBdata), pointer :: State_ptr => null()

    private nvar, na, i_de4, i_c4, amin, ai, ya, dya, Eqs_CF, ItP, GetCorrect_initial

contains

    ! Equation of state w(a) for dark energy
    real(dl) function EoS_CF(a, adotoa, grhov_t, grhoc_t)
        real(dl), intent(in) :: a, adotoa, grhoc_t, grhov_t

        select case(CFT%WForm)
        case (CPL)
            ! CPL parametrization: w(a) = w0 + w1*(1-a)
            EoS_CF = CFP%w0 + CFP%w1*(1._dl - a)
        case (HDE)
            ! Holographic Dark Energy
            EoS_CF = -1._dl/3._dl - 2*sqrt(abs(grhov_t)/3._dl)/(3*CFP%c*adotoa)
        case (NADE)
            ! New Agegraphic Dark Energy
            EoS_CF = -1._dl + 2*sqrt(abs(grhov_t)/3._dl)/(3*a*ItP%n*adotoa)
        case default
            stop 'Invalid w form in EoS_CF'
        end select
    end function EoS_CF

    ! Energy transfer coupling Q (returns 8*pi*G*a^3*Q)
    real(dl) function Coup_CF(State, a, adotoa, grhov_t, grhoc_t)
        use results, only: CAMBdata
        use classes
        class(TCAMBdata), intent(in) :: State
        real(dl), intent(in) :: a, adotoa, grhoc_t, grhov_t
        real(dl) :: wde, H0

        select case(CFT%QForm)
        case (Hrde)
            ! Q = beta * H * rho_de
            Coup_CF = CFP%beta * adotoa * grhov_t
        case (Hrc)
            ! Q = beta * H * rho_c
            Coup_CF = CFP%beta * adotoa * grhoc_t
        case (H0rde)
            ! Q = beta * H_0 * rho_de
            select type(State)
            class is (CAMBdata)
                H0 = sqrt(State%grhocrit/3._dl)
                Coup_CF = CFP%beta * H0 * a * grhov_t
            end select
        case (H0rc)
            ! Q = beta * H_0 * rho_c
            select type(State)
            class is (CAMBdata)
                H0 = sqrt(State%grhocrit/3._dl)
                Coup_CF = CFP%beta * H0 * a * grhoc_t
            end select
        case default
            ! Alternative coupling forms (e.g., NGCG model)
            wde = EoS_CF(a, adotoa, grhov_t, grhoc_t)
            Coup_CF = -3._dl * wde * CFP%beta * adotoa * grhov_t * grhoc_t / (grhov_t + grhoc_t)
        end select
    end function Coup_CF

    ! Perturbation coupling in density/pressure (returns gC = 8*pi*G*a^3*C)
    subroutine PerturCoupC_CF(State, a, adotoa, grhov_t, grhoc_t, gC)
        use classes
        class(TCAMBdata), intent(in) :: State
        real(dl), intent(in) :: a, adotoa, grhoc_t, grhov_t
        real(dl), intent(out) :: gC(3)
        real(dl) :: gQ

        gQ = Coup_CF(State, a, adotoa, grhov_t, grhoc_t)

        select case(CFT%QForm)
        case (Hrde, H0rde)
            ! Coupling proportional to rho_de
            gC(1) = gQ
            gC(2) = 0._dl
            gC(3) = 0._dl
        case (Hrc, H0rc)
            ! Coupling proportional to rho_c
            gC(1) = 0._dl
            gC(2) = gQ
            gC(3) = 0._dl
        case default
            ! NGCG or other models
            gC(1) = gQ
            gC(2) = gQ
            gC(3) = 0._dl
        end select
    end subroutine PerturCoupC_CF

    ! Perturbation coupling in velocity (returns gD = 8*pi*G*a^3*D)
    subroutine PerturCoupD_CF(State, a, adotoa, grhov_t, grhoc_t, gD)
        use classes
        class(TCAMBdata), intent(in) :: State
        real(dl), intent(in) :: a, adotoa, grhoc_t, grhov_t
        real(dl), intent(out) :: gD(2)
        real(dl) :: gQ

        gQ = Coup_CF(State, a, adotoa, grhov_t, grhoc_t)

        select case(CFT%CovQForm)
        case (uc)
            ! Coupling in dark matter rest frame
            gD(1) = 0._dl
            gD(2) = gQ
        case (ude)
            ! Coupling in dark energy rest frame
            gD(1) = gQ
            gD(2) = 0._dl
        case default
            stop 'Invalid covariant form in PerturCoupD_CF'
        end select
    end subroutine PerturCoupD_CF

    ! Determine logical flags for perturbation evolution
    subroutine LogicalPertur_CF(perturDE, evolveVc)
        logical, intent(out) :: perturDE, evolveVc
        real(dl) :: gD(2)

        ! Default: evolve perturbations, don't evolve v_c
        perturDE = .true.
        evolveVc = .false.

        ! Check if model is cosmological constant
        select case(CFT%WForm)
        case (CPL)
            if (CFP%w0 == -1._dl .and. CFP%w1 == 0._dl .and. CFP%beta == 0._dl) then
                perturDE = .false.
            end if
        case (HDE, NADE)
            perturDE = .true.
        case default
            stop 'Invalid w form in LogicalPertur_CF'
        end select

        ! Check if we need to evolve dark matter velocity
        ! Use dummy values to check coupling structure
        call PerturCoupD_CF(State_ptr, 1._dl, 1._dl, 1._dl, 1._dl, gD)
        if (gD(1) /= 0._dl) evolveVc = .true.

    end subroutine LogicalPertur_CF

    ! Set initial conditions for background evolution
    subroutine InitialCondition_CF(State, a, y)
        use results, only: CAMBdata
        use classes
        class(TCAMBdata), intent(in) :: State
        real(dl), intent(out) :: a, y(nvar)
        real(dl) :: Fa

        select type(State)
        class is (CAMBdata)
            select case(CFT%WForm)
            case (CPL, HDE)
                ! Start from today (a=1) with known densities
                a = 1._dl
                y(i_de4) = State%grhov
                y(i_c4) = State%grhoc
            case (NADE)
                ! Start from early time for NADE model
                a = amin
                Fa = 1._dl/((State%grhoc + State%grhob)/(State%grhog + State%grhornomass)*a + 1._dl)
                y(i_de4) = 0.75_dl * (State%adotrad * ItP%n * a * (1._dl + sqrt(Fa)))**2
                y(i_c4) = State%grhoc * a**(1._dl + ItP%diff)
            end select
        end select
    end subroutine InitialCondition_CF

    ! Main ODE solver: solve background equations from a=1 to a=amin
    subroutine Solve_CF(State)
        use classes
        use results, only: CAMBdata
        class(TCAMBdata), intent(in), target :: State
        real(dl) :: EV, c(24), w(nvar, 9), y(nvar), dy(nvar), tol1, a, aend
        integer :: ind, j
        external :: dverk

        ! Set module-level State pointer for use in Eqs_CF
        select type(State)
        class is (CAMBdata)
            State_ptr => State
        end select

        select case(CFT%WForm)
        case (CPL, HDE)
            a_trans = 0._dl
            ind = 1
            tol1 = 1.d-8

            call InitialCondition_CF(State, a, y)

            ! Integrate backwards in time from a=1 to a=amin
            do j = 1, na
                aend = amin**((j - 1._dl)/(na - 1._dl))
                call dverk(EV, nvar, Eqs_CF, a, y, aend, tol1, ind, c, nvar, w)
                call Eqs_CF(EV, nvar, a, y, dy)

                ! Store results (reversed order: from amin to 1)
                ai(na - j + 1) = a
                ya(na - j + 1, :) = y
                dya(na - j + 1, :) = dy

                ! Track transition point where rho_de goes negative
                if (y(i_de4) <= 0._dl .and. a_trans == 0._dl) a_trans = a
            end do

        case (NADE)
            a_trans = 0._dl
            ind = 1
            tol1 = 1.d-8

            ! Find correct initial conditions via iteration
            call GetCorrect_initial(State)
            call InitialCondition_CF(State, a, y)

            ! Integrate forward in time from a=amin to a=1
            do j = 1, na
                aend = amin**((na - j)/(na - 1._dl))
                call dverk(EV, nvar, Eqs_CF, a, y, aend, tol1, ind, c, nvar, w)
                call Eqs_CF(EV, nvar, a, y, dy)

                ai(j) = a
                ya(j, :) = y
                dya(j, :) = dy
            end do
        end select

    end subroutine Solve_CF

    ! Background evolution equations: dy/da
    subroutine Eqs_CF(State_EV, n, a, y, dy)
        use results, only: CAMBdata
        integer, intent(in) :: n
        real(dl), intent(in) :: State_EV, a, y(n)
        real(dl), intent(out) :: dy(n)
        real(dl) :: a2, adotoa, grhov_t, grhoc_t, gQ, wde

        a2 = a*a
        grhov_t = y(i_de4)/a2
        grhoc_t = y(i_c4)/a2

        adotoa = ConformalH(State_ptr, a, grhov_t, grhoc_t)
        wde = EoS_CF(a, adotoa, grhov_t, grhoc_t)
        gQ = Coup_CF(State_ptr, a, adotoa, grhov_t, grhoc_t)

        ! Evolution equations:
        ! d(rho_de*a^2)/da = (1 - 3*w)*rho_de*a^2/a + Q*a^2/(H*a)
        dy(i_de4) = (1 - 3*wde)*grhov_t*a + gQ*a/adotoa

        ! d(rho_c*a^2)/da = rho_c*a^2/a - Q*a^2/(H*a)
        dy(i_c4) = grhoc_t*a - gQ*a/adotoa

    end subroutine Eqs_CF

    ! Get interpolated background densities at scale factor a
    subroutine Get_grhodea2_grhoca2(a, grhov_t, grhoc_t)
        real(dl), intent(in) :: a
        real(dl), intent(out) :: grhov_t, grhoc_t
        integer :: klo

        ! Find interpolation index
        klo = int(na - (na - 1)*log(a)/log(amin))
        klo = max(min(klo, na - 1), 1)

        ! Hermite interpolation
        grhov_t = Hermite(ai, ya(:, i_de4), dya(:, i_de4), a, klo)/a**2
        grhoc_t = Hermite(ai, ya(:, i_c4), dya(:, i_c4), a, klo)/a**2

    end subroutine Get_grhodea2_grhoca2

    ! Find correct initial conditions for NADE model via iteration
    subroutine GetCorrect_initial(State)
        use classes
        use results, only: CAMBdata
        class(TCAMBdata), intent(in), target :: State
        integer, parameter :: n = 2, itmax = 100
        real(dl), dimension(n) :: x0, x1, f0, f1, dx, df, v1, v2
        real(dl), dimension(n, n) :: D0, H0, H1
        real(dl) :: t1, t2, tol1, delta
        integer :: i

        tol1 = 1.d-5
        delta = 1.d-3
        x0(1) = ItP%n
        x0(2) = ItP%diff

        ! Get initial derivative matrix
        call IteraFcns(State, x0, f0)

        x1(1) = x0(1) + delta
        x1(2) = x0(2)
        call IteraFcns(State, x1, f1)
        D0(:, 1) = (f1 - f0)/delta

        x1(1) = x0(1)
        x1(2) = x0(2) + delta
        call IteraFcns(State, x1, f1)
        D0(:, 2) = (f1 - f0)/delta

        ! Invert derivative matrix
        H0(1, 1) = D0(2, 2)
        H0(2, 1) = -D0(2, 1)
        H0(1, 2) = -D0(1, 2)
        H0(2, 2) = D0(1, 1)
        H0 = H0/(D0(1, 1)*D0(2, 2) - D0(1, 2)*D0(2, 1))

        ! Broyden iteration
        do i = 1, itmax
            if (sqrt(sum(f0**2)) < tol1) exit

            x1 = x0 - matmul(H0, f0)
            call IteraFcns(State, x1, f1)

            dx = x1 - x0
            df = f1 - f0
            v1 = matmul(H0, df)
            v2 = dx - v1
            t1 = dot_product(v2, dx)
            t2 = dot_product(dx, v1)
            H1 = (1 + t1/t2)*H0

            x0 = x1
            f0 = f1
            H0 = H1

            if (i >= itmax) stop 'BroydenIteration not converge in GetCorrect_initial'
        end do

    contains

        subroutine IteraFcns(State, x, f)
            use classes
            use results, only: CAMBdata
            class(TCAMBdata), intent(in) :: State
            real(dl), dimension(:), intent(in) :: x
            real(dl), dimension(:), intent(out) :: f
            real(dl) :: EV, c(24), w(nvar, 9), y(nvar), tol1, a
            integer :: ind
            external :: dverk

            ind = 1
            tol1 = 1.d-8
            ItP%n = x(1)
            ItP%diff = x(2)

            call InitialCondition_CF(State, a, y)
            call dverk(EV, nvar, Eqs_CF, a, y, a, tol1, ind, c, nvar, w)  ! Handle first point
            call dverk(EV, nvar, Eqs_CF, a, y, 1._dl, tol1, ind, c, nvar, w)  ! Evolve to a=1

            ! Residuals: difference from target values at a=1
            select type(State)
            class is (CAMBdata)
                f(1) = (y(i_de4) - State%grhov)/State%grhov
                f(2) = (y(i_c4) - State%grhoc)/State%grhoc
            end select
        end subroutine IteraFcns

    end subroutine GetCorrect_initial

end module CoupledFluidModels
