! ===========================================================================================================
! ===========================================================================================================

MODULE donnees

    use math

    implicit none

    real(rp), parameter :: gamma_eau_air = 0.073_rp

contains

    function pression(x, centreBulle, rayonBulle)
        ! paramètres
        real(rp), dimension(2), intent(in) :: x, centreBulle
        real(rp), intent(in) :: rayonBulle

        ! return
        real(rp) :: pression

        ! variables locales
        real(rp) :: t
        real(rp) :: max_P

        t = fnorme2(x - centreBulle)

        if (t < rayonBulle) then
            max_P = 2.0_rp * gamma_eau_air / rayonBulle
            pression = (-max_P / rayonBulle) * t + max_P
        else
            pression = 0.0_rp
        end if
    end function

    subroutine init_pression(x, centreBulle, rayonBulle, P)
        ! paramètres
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: centreBulle
        real(rp), intent(in) :: rayonBulle
        real(rp), dimension(:), intent(out) :: P

        ! variables locales
        integer :: np, i

        np = size(P)
        do i = 1, np
            P(i) = pression(x(i, :), centreBulle, rayonBulle)
        end do
    end subroutine

END MODULE donnees
