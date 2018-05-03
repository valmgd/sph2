! ===========================================================================================================
! ===========================================================================================================

MODULE donnees

    use math
    use sph

    implicit none

    real(rp), parameter :: gamma_eau_air = 0.073_rp

contains

    ! -------------------------------------------------------------------------------------------------------
    ! lire valeur en cherchant balises
    ! -------------------------------------------------------------------------------------------------------
    subroutine readValues(nom_fichier, n, bornes, centre, rayon)
        ! paramètres
        character(len=*), intent(in) :: nom_fichier
        integer, intent(out) :: n
        real(rp), dimension(:, :), intent(out) :: bornes
        real(rp), dimension(:), intent(out) :: centre
        real(rp), intent(out) :: rayon

        ! variables locales
        character(len=20) :: ligne
        integer :: i

        ligne = " "
        open(unit = 10, file = nom_fichier)
        do while (trim(ligne) /= "#n")
            read (10, *) ligne
        end do
        read (10, *) n
        close(10)

        ligne = " "
        open(unit = 10, file = nom_fichier)
        do while (trim(ligne) /= "#bornes")
            read (10, *) ligne
        end do
        do i = 1, d
            read (10, *) bornes(i, :)
        end do
        close(10)

        ligne = " "
        open(unit = 10, file = nom_fichier)
        do while (trim(ligne) /= "#centre")
            read (10, *) ligne
        end do
        read (10, *) centre
        close(10)

        ligne = " "
        open(unit = 10, file = nom_fichier)
        do while (trim(ligne) /= "#rayon")
            read (10, *) ligne
        end do
        read (10, *) rayon
        close(10)
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! pression, fonction affine, nulle au bord, ayant une valeur max au centre de la bulle
    ! -------------------------------------------------------------------------------------------------------
    ! x : point dans la bulle auquel on veut estimer la pression
    ! centreBulle : centre
    ! rayonBulle : rayon
    function pression(x, centreBulle, rayonBulle)
        ! paramètres
        real(rp), dimension(d), intent(in) :: x, centreBulle
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



    ! -------------------------------------------------------------------------------------------------------
    ! initialise la pression pour toutes les particules
    ! -------------------------------------------------------------------------------------------------------
    ! part : liste de particules
    subroutine init_pression_bulle(centreBulle, rayonBulle, part)
        ! paramètres
        real(rp), dimension(:), intent(in) :: centreBulle
        real(rp), intent(in) :: rayonBulle
        type(Particules), intent(inout) :: part

        ! variables locales
        integer :: i

        do i = 1, part%n
            part%P(i) = pression(part%x(i, :), centreBulle, rayonBulle)
        end do
    end subroutine

    subroutine init_pression_pave(bornes, part)
        ! paramètres
        real(rp), dimension(:, :), intent(in) :: bornes
        type(Particules), intent(inout) :: part

        ! variables locales

        part%P = 0.0_rp
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! maille un carré / cube / hypercube... et écrit les valeurs dans un fichier
    ! -------------------------------------------------------------------------------------------------------
    ! d_Omega : dimension du pavé à mailler
    ! n : nombre de points de maillage par axe du domaine (n points en x, n points en y etc)
    ! bornes : bornes du pavé à mailler -> xmin, xmax \n ymin, ymax \n zmin, zmax
    ! nom_fichier : début du nom des fichiers dans lequel sont écrits les points et l'enveloppe du pavé
    ! part : liste de particules retournées (coordonnées + volumes)
    subroutine pave(d_Omega, n, bornes, nom_fichier, part)
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
        allocate(part%w(part%n), part%rho(part%n), part%u(part%n, d_Omega), part%P(part%n))
        part%w = dx**d_Omega
        part%R = I_SPH * dx

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
        allocate(part%w(part%n), part%rho(part%n), part%u(part%n, d_Omega), part%P(part%n))
        part%w = dx**d_Omega
        part%R = I_SPH * dx

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



    ! -------------------------------------------------------------------------------------------------------
    ! calcul GR(1)
    ! -------------------------------------------------------------------------------------------------------
    subroutine normale_surface_GR(part, nom_fichier)
        ! paramètres
        type(Particules), intent(in) :: part
        character(len=*), intent(in) :: nom_fichier

        ! variables locales
        integer :: i, k
        real(rp), dimension(d) :: grad
        real(rp), dimension(part%n, d * 2) :: nvec
        real(rp) :: length

        length = 0.5_rp * (fnorme2(part%x(part%n/2, :) - part%x(part%n/2+1, :)))

        do i = 1, part%n
            call GR(i, part, (/ (1.0_rp, k=1, part%n) /), grad)
            if (fnorme2(grad) > 10.0**(-10)) then
                nvec(i, :) = (/ part%x(i, :), (grad/fnorme2(grad))*length /)
            else
                nvec(i, :) = 0.0_rp
            end if
        end do

        call writeMat(nvec, nom_fichier)
    end subroutine

END MODULE donnees
