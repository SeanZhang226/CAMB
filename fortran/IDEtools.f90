! IDEtools module for Interacting Dark Energy
! Contains utility functions for IDE calculations
! Ported from liaocrane/IDECAMB

module IDEtools
    use precision
    use MassiveNu
    implicit none

contains

    ! Calculate conformal Hubble parameter H*a at scale factor a
    ! adotoa = H*a = sqrt((8*pi*G/3) * rho_total)
    function ConformalH(State, a, grhov_t, grhoc_t) result(adotoa)
        use classes
        use results, only: CAMBdata
        class(TCAMBdata), intent(in), target :: State
        real(dl), intent(in) :: a, grhov_t, grhoc_t
        real(dl) :: adotoa, a2, grho, rhonu
        integer :: nu_i

        select type (State)
        class is (CAMBdata)
            a2 = a*a
            ! Total energy density: curvature + baryons + dark energy + dark matter + photons + massless neutrinos
            grho = State%grhok + State%grhob/a + grhov_t + grhoc_t + (State%grhog + State%grhornomass)/a2

            ! Add massive neutrinos if present
            if (State%CP%Num_Nu_massive /= 0) then
                do nu_i = 1, State%CP%Nu_mass_eigenstates
                    call ThermalNuBackground%rho(a*State%nu_masses(nu_i), rhonu)
                    grho = grho + rhonu*State%grhormass(nu_i)/a2
                end do
            end if

            ! H*a = sqrt(8*pi*G*rho/3)
            adotoa = sqrt(grho/3._dl)
        end select
    end function ConformalH

    ! Hermite interpolation using function values ya and derivatives dya
    ! Interpolates between xa(klo) and xa(klo+1) to find y(x)
    function Hermite(xa, ya, dya, x, klo) result(y)
        real(dl), dimension(:), intent(in) :: xa, ya, dya
        real(dl), intent(in) :: x
        integer, intent(in) :: klo
        real(dl) :: y, alpha, beta, h
        integer :: khi

        khi = klo + 1
        h = xa(khi) - xa(klo)
        alpha = (xa(khi) - x)/h
        beta = (x - xa(klo))/h

        ! Hermite cubic interpolation formula
        y = ya(klo)*(1 + 2*beta)*alpha**2 + ya(khi)*(1 + 2*alpha)*beta**2 + &
            dya(klo)*beta*h*alpha**2 - dya(khi)*alpha*h*beta**2

    end function Hermite

end module IDEtools
