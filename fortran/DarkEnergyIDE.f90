! Interacting Dark Energy (IDE) Model Module
! Based on IDECAMB (liaocrane/IDECAMB)
! Simplified implementation for coupled dark energy-dark matter interaction

module DarkEnergyIDE
    use precision
    use DarkEnergyInterface
    use classes
    implicit none
    private

    type, extends(TDarkEnergyEqnOfState) :: TDarkEnergyIDE
        ! Coupling parameters
        real(dl) :: xi_de = 0._dl  ! Coupling constant proportional to dark energy density
        real(dl) :: xi_c = 0._dl   ! Coupling constant proportional to dark matter density
        
        ! Background evolution storage
        real(dl), allocatable :: a_table(:), grhov_table(:), grhoc_table(:)
        integer :: n_table
        logical :: tables_initialized = .false.
        
    contains
        procedure :: ReadParams => TDarkEnergyIDE_ReadParams
        procedure :: Init => TDarkEnergyIDE_Init
        procedure :: PrintFeedback => TDarkEnergyIDE_PrintFeedback
        procedure :: BackgroundDensityAndPressure => TDarkEnergyIDE_BackgroundDensityAndPressure
        procedure :: PerturbedStressEnergy => TDarkEnergyIDE_PerturbedStressEnergy
        procedure :: diff_rhopi_Add_Term => TDarkEnergyIDE_diff_rhopi_Add_Term
        procedure :: PerturbationEvolve => TDarkEnergyIDE_PerturbationEvolve
        procedure :: PerturbationInitial => TDarkEnergyIDE_PerturbationInitial
    end type TDarkEnergyIDE

    public TDarkEnergyIDE
    
contains

    subroutine TDarkEnergyIDE_ReadParams(this, Ini)
        use IniObjects
        class(TDarkEnergyIDE) :: this
        class(TIniFile), intent(in) :: Ini
        
        ! Read standard w, wa parameters first
        call this%TDarkEnergyEqnOfState%ReadParams(Ini)
        
        ! Read IDE coupling parameters
        this%xi_de = Ini%Read_Double('xi_de', 0._dl)
        this%xi_c = Ini%Read_Double('xi_c', 0._dl)
        
    end subroutine TDarkEnergyIDE_ReadParams

    subroutine TDarkEnergyIDE_Init(this, State)
        class(TDarkEnergyIDE), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State
        
        ! Initialize parent class
        call this%TDarkEnergyEqnOfState%Init(State)
        
        ! Check if we have coupling (non-zero xi_de or xi_c)
        if (abs(this%xi_de) > 1e-10 .or. abs(this%xi_c) > 1e-10) then
            this%is_cosmological_constant = .false.
            this%num_perturb_equations = 2  ! clxq and vq for dark energy
        else
            ! No coupling, behave like standard dark energy
            this%is_cosmological_constant = .not. this%use_tabulated_w .and. &
                abs(this%w_lam + 1._dl) < 1.e-6_dl .and. this%wa==0._dl
            this%num_perturb_equations = 0
        end if
        
    end subroutine TDarkEnergyIDE_Init

    subroutine TDarkEnergyIDE_PrintFeedback(this, FeedbackLevel)
        class(TDarkEnergyIDE) :: this
        integer, intent(in) :: FeedbackLevel
        
        call this%TDarkEnergyEqnOfState%PrintFeedback(FeedbackLevel)
        
        if (FeedbackLevel > 0) then
            write(*,'("IDE coupling: (xi_de, xi_c) = (", f10.6,", ", f10.6, ")")') &
                this%xi_de, this%xi_c
        end if
        
    end subroutine TDarkEnergyIDE_PrintFeedback

    subroutine TDarkEnergyIDE_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
        class(TDarkEnergyIDE), intent(inout) :: this
        real(dl), intent(in) :: grhov, a
        real(dl), intent(out) :: grhov_t
        real(dl), optional, intent(out) :: w
        real(dl) :: coupling_factor, adotoa_est
        
        ! For small couplings, use perturbative approach
        if (abs(this%xi_de) < 1e-3 .and. abs(this%xi_c) < 1e-3) then
            ! Use standard evolution as base
            call this%TDarkEnergyEqnOfState%BackgroundDensityAndPressure(grhov, a, grhov_t, w)
            
            ! Apply small coupling correction
            if (abs(this%xi_de) > 1e-10 .or. abs(this%xi_c) > 1e-10) then
                ! Approximate correction: Q ~ xi * H * rho
                ! This modifies the density slightly
                adotoa_est = sqrt(grhov_t / (3._dl * a**2))
                coupling_factor = 1._dl - (this%xi_de - this%xi_c) * log(a)
                grhov_t = grhov_t * coupling_factor
            end if
        else
            ! For larger couplings, use standard evolution
            ! (More sophisticated integration would go here in production)
            call this%TDarkEnergyEqnOfState%BackgroundDensityAndPressure(grhov, a, grhov_t, w)
        end if
        
    end subroutine TDarkEnergyIDE_BackgroundDensityAndPressure

    subroutine TDarkEnergyIDE_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        class(TDarkEnergyIDE), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix
        real(dl) :: clxq, vq, coupling_term
        
        ! Initialize to zero
        dgrhoe = 0._dl
        dgqe = 0._dl
        
        ! If no coupling, return
        if (abs(this%xi_de) < 1e-10 .and. abs(this%xi_c) < 1e-10) return
        
        ! If we have perturbation equations
        if (this%num_perturb_equations > 0 .and. w_ix > 0) then
            clxq = ay(w_ix)      ! Dark energy density perturbation
            vq = ay(w_ix + 1)    ! Dark energy velocity perturbation
            
            ! Dark energy perturbations
            dgrhoe = grhov_t * clxq
            dgqe = grhov_t * vq
            
            ! Add coupling term to dark matter perturbations
            ! This represents the interaction Q
            coupling_term = (this%xi_de * grhov_t + this%xi_c * dgrho) * adotoa / k
            dgqe = dgqe + coupling_term
        end if
        
    end subroutine TDarkEnergyIDE_PerturbedStressEnergy

    function TDarkEnergyIDE_diff_rhopi_Add_Term(this, dgrhoe, dgqe, grho, gpres, w, grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
        class(TDarkEnergyIDE), intent(in) :: this
        real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, grhok, w, adotoa, &
            k, grhov_t, z, k2, yprime(:), y(:), Kf1
        integer, intent(in) :: w_ix
        real(dl) :: ppiedot
        
        ppiedot = 0._dl
        
        ! Add anisotropic stress contribution if needed
        ! For simple IDE model, we assume no anisotropic stress
        
    end function TDarkEnergyIDE_diff_rhopi_Add_Term

    subroutine TDarkEnergyIDE_PerturbationEvolve(this, ayprime, w, w_ix, a, adotoa, k, z, y)
        class(TDarkEnergyIDE), intent(in) :: this
        real(dl), intent(inout) :: ayprime(:)
        real(dl), intent(in) :: a, adotoa, k, z, y(:), w
        integer, intent(in) :: w_ix
        real(dl) :: clxq, vq, clxc, vc, grhoc_frac
        real(dl) :: coupling_de, coupling_c, cs2_eff
        
        if (this%num_perturb_equations == 0 .or. w_ix == 0) return
        
        ! Get perturbation variables
        clxq = y(w_ix)       ! DE density perturbation
        vq = y(w_ix + 1)     ! DE velocity perturbation
        
        ! Effective sound speed
        cs2_eff = this%cs2_lam
        
        ! Evolution equations for dark energy perturbations with coupling
        ! d(clxq)/dtau = -3 * adotoa * (cs2_eff - w) * clxq - (1 + w) * k * vq + coupling
        coupling_de = this%xi_de * adotoa * 3._dl * (1._dl + w)
        ayprime(w_ix) = -3._dl * adotoa * (cs2_eff - w) * clxq - (1._dl + w) * k * vq + coupling_de * clxq
        
        ! d(vq)/dtau = -adotoa * (1 - 3*cs2_eff) * vq + k * cs2_eff * clxq / (1 + w) + coupling
        coupling_c = this%xi_c * adotoa
        ayprime(w_ix + 1) = -adotoa * (1._dl - 3._dl * cs2_eff) * vq + &
            k * cs2_eff * clxq / (1._dl + w) + coupling_c * vq
        
    end subroutine TDarkEnergyIDE_PerturbationEvolve

    subroutine TDarkEnergyIDE_PerturbationInitial(this, y, a, tau, k)
        class(TDarkEnergyIDE), intent(in) :: this
        real(dl), intent(out) :: y(:)
        real(dl), intent(in) :: a, tau, k
        
        ! Initialize perturbations to zero (adiabatic initial conditions)
        y = 0._dl
        
    end subroutine TDarkEnergyIDE_PerturbationInitial

end module DarkEnergyIDE
