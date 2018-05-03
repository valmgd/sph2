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

    ! -------------------------------------------------------------------------------------------------------
    ! maille un carré / cube / hypercube... et écrit les valeurs dans un fichier
    ! -------------------------------------------------------------------------------------------------------
    ! d_Omega : dimension du pavé à mailler
    ! n : nombre de points de maillage par axe du domaine (n points en x, n points en y etc)
    ! bornes : bornes du pavé à mailler -> xmin, xmax \n ymin, ymax \n zmin, zmax
    ! nom_fichier : début du nom des fichiers dans lequel sont écrits les points et l'enveloppe du pavé
    ! part : liste de particules retournées (coordonnées + volumes)
    subroutine cube(d_Omega, n, bornes, nom_fichier, part)
        ! paramètres
        integer, intent(in) :: d_Omega, n
        real(rp), dimension(d_Omega, 2), intent(in) :: bornes
        character(len=*), intent(in) :: nom_fichier
        type(Particules), intent(out) :: part

        ! variables locales
        real(rp), dimension(:, :), allocatable :: subd_axes
        integer :: i, np
        real(rp) :: dx

        dx = (bornes(1, 2) - bornes(1, 1)) / n
        allocate(subd_axes(n, d_Omega))
        do i = 1, d_Omega
            subd_axes(:, i) = linspace(bornes(i, 1) + dx/2.0_rp, bornes(i, 2) - dx/2.0_rp, n)
        end do

        call meshgrid(subd_axes, part%x)
        deallocate(subd_axes)

        part%n = size(part%x, 1)
        allocate(part%w(part%n))
        part%w = dx**d_Omega

        call writeMat(part%x, nom_fichier // "_points.dat")
        open(unit = 10, file = nom_fichier // "_d.dat")
        write (10, *) d_Omega
        close(10)
        if (d_Omega == 2) then
            open(unit = 10, file = nom_fichier // "_enveloppe.dat")
            write (10, *) bornes(1, 1), bornes(2, 1)
            write (10, *) bornes(1, 2), bornes(2, 1)
            write (10, *) bornes(1, 2), bornes(2, 2)
            write (10, *) bornes(1, 1), bornes(2, 2)
            write (10, *) bornes(1, 1), bornes(2, 1)
            close(10)
        else if (d_Omega == 3) then
            open(unit = 10, file = nom_fichier // "_enveloppe.dat")
            write (10, *) bornes(1, 1), bornes(2, 1), bornes(3, 1)
            write (10, *) bornes(1, 1), bornes(2, 1), bornes(3, 2)
            write (10, *) bornes(1, 1), bornes(2, 2), bornes(3, 2)
            write (10, *) bornes(1, 1), bornes(2, 2), bornes(3, 1)
            write (10, *) bornes(1, 1), bornes(2, 1), bornes(3, 1)
            write (10, *)

            write (10, *) bornes(1, 2), bornes(2, 1), bornes(3, 1)
            write (10, *) bornes(1, 2), bornes(2, 1), bornes(3, 2)
            write (10, *) bornes(1, 2), bornes(2, 2), bornes(3, 2)
            write (10, *) bornes(1, 2), bornes(2, 2), bornes(3, 1)
            write (10, *) bornes(1, 2), bornes(2, 1), bornes(3, 1)
            close(10)
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! maille une bulle et écrit les valeurs dans un fichier
    ! -------------------------------------------------------------------------------------------------------
    ! d_Omega : dimension du pavé à mailler
    ! n : nombre de points de maillage par axe du domaine (n points en x, n points en y etc)
    ! bornes : bornes du pavé à mailler -> xmin, xmax \n ymin, ymax \n zmin, zmax
    ! nom_fichier : début du nom des fichiers dans lequel sont écrits les points et l'enveloppe du pavé
    ! part : liste de particules retournées (coordonnées + volumes)
    subroutine bulle(d_Omega, n, centre, rayon, nom_fichier, part)
        ! paramètres
        integer, intent(in) :: d_Omega, n
        real(rp), dimension(d_Omega), intent(in) :: centre
        real(rp), intent(in) :: rayon
        character(len=*), intent(in) :: nom_fichier
        type(Particules), intent(out) :: part

        ! variables locales
        real(rp), dimension(:, :), allocatable :: subd_axes
        integer :: i, np
        real(rp) :: dx
        real(rp), dimension(100) :: t
        real(rp), dimension(100, d_Omega) :: cercle

        call meshCircle(d_Omega, centre, rayon, n, part%x)
        dx = (part%x(2, 1) - part%x(1, 1))
        part%n = size(part%x, 1)
        allocate(part%w(part%n))
        part%w = dx**d_Omega

        ! sauvegarde particules et dimension
        call writeMat(part%x, nom_fichier // "_points.dat")
        open(unit = 10, file = nom_fichier // "_d.dat")
        write (10, *) d_Omega
        close(10)

        ! construction enveloppe
        t = linspace(0.0_rp, 2.0_rp * pi, 100)
        do i = 1, 100
            cercle(i, :) = centre + rayon * (/ cos(t(i)), sin(t(i)) /)
        end do

        ! sauvegarde enveloppe
        if (d_Omega == 2) then
            call writeMat(cercle, nom_fichier // "_enveloppe.dat")
        else if (d_Omega == 3) then
        open(unit = 10, file = nom_fichier // "_enveloppe.dat")
        write (10, *)
        close(10)
        end if
    end subroutine

END MODULE donnees
