! ===========================================================================================================
! ===========================================================================================================

MODULE donnees

    use math
    use sph

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

    subroutine cube(d_Omega, n, xmin, xmax, nom_fichier, x, w)
        ! paramètres
        integer, intent(in) :: d_Omega, n
        real(rp), intent(in) :: xmax, xmin
        character(len=*), intent(in) :: nom_fichier
        real(rp), dimension(:, :), allocatable, intent(out) :: x
        real(rp), dimension(:), allocatable, intent(out) :: w

        ! variables locales
        real(rp), dimension(:, :), allocatable :: subd_axes
        integer :: i, np
        real(rp) :: dx

        dx = (xmax - xmin) / n
        allocate(subd_axes(n, d_Omega))
        do i = 1, d_Omega
            subd_axes(:, i) = linspace(xmin + dx/2.0_rp, xmax - dx/2.0_rp, n)
        end do

        call meshgrid(subd_axes, x)
        deallocate(subd_axes)

        np = size(x, 1)
        allocate(w(np))
        w = dx

        call writeMat(x, nom_fichier // "_points.dat")
        open(unit = 10, file = nom_fichier // "_enveloppe.dat")
        write (10, *) xmin, xmin
        write (10, *) xmax, xmin
        write (10, *) xmax, xmax
        write (10, *) xmin, xmax
        write (10, *) xmin, xmin
        close(10)
    end subroutine

END MODULE donnees
