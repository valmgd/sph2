! ===========================================================================================================
! Module contenant les fonctions suivantes :
!
!   subroutine set_DONNEES_SIGMA(valeur)
!   subroutine readValues(nom_fichier, d_Omega, sigma, intervalle, n, bornes, centre, rayon)
!   subroutine pave(d_Omega, n, bornes, nom_fichier, part)
!   subroutine bulle(d_Omega, n, centre, rayon, nom_fichier, part)
!   subroutine init_var_bulle(centreBulle, rayonBulle, part)
!   subroutine init_var_pave(bornes, part)
!   subroutine write_var(part, f_rho, f_u, f_P)
!   subroutine normale_surface_GR(part, nom_fichier)
!   subroutine quarter(part, centre)
! ===========================================================================================================

MODULE donnees

    ! libraries
    use sph

    implicit none

    real(rp), save :: DONNEES_SIGMA

contains

    subroutine set_DONNEES_SIGMA(valeur)
        ! paramètres
        real(rp), intent(in) :: valeur
        DONNEES_SIGMA = valeur
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! lire valeur en cherchant balises
    ! -------------------------------------------------------------------------------------------------------
    subroutine readValues(nom_fichier, d_Omega, sigma, intervalle, n, bornes_ext, bornes_int, centre, rayon)
        ! paramètres
        character(len=*), intent(in) :: nom_fichier
        integer, intent(out) :: d_Omega
        real(rp), intent(out) :: sigma
        real(rp), intent(out) :: intervalle
        integer, intent(out) :: n
        real(rp), dimension(:, :), allocatable, intent(out) :: bornes_ext, bornes_int
        real(rp), dimension(:), allocatable, intent(out) :: centre
        real(rp), intent(out) :: rayon

        ! variables locales
        character(len=20) :: ligne
        integer :: i

        ligne = " "
        open(unit = 10, file = nom_fichier)
        do while (trim(ligne) /= "#dimension")
            read (10, *) ligne
        end do
        read (10, *) d_Omega
        close(10)

        ligne = " "
        open(unit = 10, file = nom_fichier)
        do while (trim(ligne) /= "#sigma")
            read (10, *) ligne
        end do
        read (10, *) sigma
        close(10)

        ligne = " "
        open(unit = 10, file = nom_fichier)
        do while (trim(ligne) /= "#intervalle")
            read (10, *) ligne
        end do
        read (10, *) intervalle
        close(10)

        ! affectation constantes
        call set_SPH_D(d_Omega)
        call set_DONNEES_SIGMA(sigma)
        call set_SPH_I(intervalle)
        allocate(bornes_ext(SPH_D, 2), bornes_int(SPH_D, 2), centre(SPH_D))

        ligne = " "
        open(unit = 10, file = nom_fichier)
        do while (trim(ligne) /= "#n")
            read (10, *) ligne
        end do
        read (10, *) n
        close(10)

        ligne = " "
        open(unit = 10, file = nom_fichier)
        do while (trim(ligne) /= "#bornes_ext")
            read (10, *) ligne
        end do
        do i = 1, SPH_D
            read (10, *) bornes_ext(i, :)
        end do
        close(10)

        ligne = " "
        open(unit = 10, file = nom_fichier)
        do while (trim(ligne) /= "#bornes_int")
            read (10, *) ligne
        end do
        do i = 1, SPH_D
            read (10, *) bornes_int(i, :)
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
    ! Check si une particule est dans un carré délimité par ses bornes.
    ! -------------------------------------------------------------------------------------------------------
    function set_color(x, bornes) result(couleur)
        ! paramètres
        real(rp), dimension(SPH_D), intent(in) :: x
        real(rp), dimension(SPH_D, 2), intent(in) :: bornes

        ! return
        integer :: couleur

        ! variables locales
        integer :: i

        couleur = 1

        do i = 1, SPH_D
            if (.not. ((bornes(i, 1) < x(i)) .and. (x(i) < bornes(i, 2)))) then
                couleur = 0
                exit
            end if
        end do
    end function

    function isInSquare(x, bornes) result(answer)
        ! paramètres
        real(rp), dimension(SPH_D), intent(in) :: x
        real(rp), dimension(SPH_D, 2), intent(in) :: bornes

        ! return
        logical :: answer

        ! variables locales
        integer :: i

        answer = .TRUE.

        do i = 1, SPH_D
            if (.not. ((bornes(i, 1) < x(i)) .and. (x(i) < bornes(i, 2)))) then
                answer = .FALSE.
                exit
            end if
        end do
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! maille un carré / cube / hypercube... et écrit les valeurs dans un fichier
    ! -------------------------------------------------------------------------------------------------------
    ! d_Omega : dimension du pavé à mailler
    ! n : nombre de points de maillage par axe du domaine (n points en x, n points en y etc)
    ! bornes : bornes du pavé à mailler -> xmin, xmax \n ymin, ymax \n zmin, zmax
    ! nom_fichier : début du nom des fichiers dans lequel sont écrits les points et l'enveloppe du pavé
    ! part : liste de particules retournées (coordonnées + volumes)
    subroutine pave(d_Omega, n, bornes_ext, bornes_int, nom_fichier, part)
        ! paramètres
        integer, intent(in) :: d_Omega, n
        real(rp), dimension(d_Omega, 2), intent(in) :: bornes_ext, bornes_int
        character(len=*), intent(in) :: nom_fichier
        type(Particules), intent(out) :: part

        ! variables locales
        real(rp), dimension(:, :), allocatable :: subd_axes
        integer :: i, np
        real(rp) :: dx

        character(len=20) :: ligne

        dx = (bornes_ext(1, 2) - bornes_ext(1, 1)) / n
        part%dx = dx
        allocate(subd_axes(n, d_Omega))
        do i = 1, d_Omega
            subd_axes(:, i) = linspace(bornes_ext(i, 1) + dx/2.0_rp, bornes_ext(i, 2) - dx/2.0_rp, n)
        end do

        ! meshgrid s'occupe de l'allocation de part%x
        call meshgrid(subd_axes, part%x)
        deallocate(subd_axes)

        part%n = size(part%x, 1)
        allocate(part%w(part%n), part%rho(part%n), part%u(part%n, d_Omega), part%P(part%n))
        allocate(part%gradR(part%n, d_Omega), part%fts(part%n, d_Omega), part%dWij(part%n, part%n, SPH_D))
        allocate(part%div_nor(part%n), part%nor(part%n, d_Omega), part%centre(d_Omega), part%color(part%n))
        part%w = dx**d_Omega
        part%R = SPH_I * dx


        open(unit = 10, file = nom_fichier // "_d.dat")
        write (10, *) d_Omega
        close(10)
        if (d_Omega == 2) then
            open(unit = 10, file = nom_fichier // "_enveloppe.dat")
            write (10, *) bornes_ext(1, 1), bornes_ext(2, 1)
            write (10, *) bornes_ext(1, 2), bornes_ext(2, 1)
            write (10, *) bornes_ext(1, 2), bornes_ext(2, 2)
            write (10, *) bornes_ext(1, 1), bornes_ext(2, 2)
            write (10, *) bornes_ext(1, 1), bornes_ext(2, 1)
            close(10)
            open(unit = 10, file = nom_fichier // "_enveloppe2.dat")
            write (10, *) bornes_int(1, 1), bornes_int(2, 1)
            write (10, *) bornes_int(1, 2), bornes_int(2, 1)
            write (10, *) bornes_int(1, 2), bornes_int(2, 2)
            write (10, *) bornes_int(1, 1), bornes_int(2, 2)
            write (10, *) bornes_int(1, 1), bornes_int(2, 1)
            close(10)
        else if (d_Omega == 3) then
            open(unit = 10, file = nom_fichier // "_enveloppe.dat")
            write (10, *) bornes_ext(1, 1), bornes_ext(2, 1), bornes_ext(3, 1)
            write (10, *) bornes_ext(1, 1), bornes_ext(2, 1), bornes_ext(3, 2)
            write (10, *) bornes_ext(1, 1), bornes_ext(2, 2), bornes_ext(3, 2)
            write (10, *) bornes_ext(1, 1), bornes_ext(2, 2), bornes_ext(3, 1)
            write (10, *) bornes_ext(1, 1), bornes_ext(2, 1), bornes_ext(3, 1)
            write (10, *)

            write (10, *) bornes_ext(1, 2), bornes_ext(2, 1), bornes_ext(3, 1)
            write (10, *) bornes_ext(1, 2), bornes_ext(2, 1), bornes_ext(3, 2)
            write (10, *) bornes_ext(1, 2), bornes_ext(2, 2), bornes_ext(3, 2)
            write (10, *) bornes_ext(1, 2), bornes_ext(2, 2), bornes_ext(3, 1)
            write (10, *) bornes_ext(1, 2), bornes_ext(2, 1), bornes_ext(3, 1)
            close(10)
        end if

        ! affectation de la couleur des particules (intérieur 0 extérieur 1)
        do i = 1, part%n
            part%color(i) = set_color(part%x(i, :), bornes_int)
        end do

        open(unit = 10, file = nom_fichier // "_points.dat")
        do i = 1, part%n
            write (10, *), part%x(i, :), part%color(i)
        end do
        close(10)

        call writeMat(transpose(bornes_ext), nom_fichier // "_scale.dat")
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
        !dx = (part%x(2, 1) - part%x(1, 1))
        dx = 2.0_rp * rayon / real(n, rp)
        part%dx = dx
        part%n = size(part%x, 1)
        allocate(part%w(part%n), part%rho(part%n), part%u(part%n, d_Omega), part%P(part%n))
        allocate(part%gradR(part%n, d_Omega), part%fts(part%n, d_Omega), part%dWij(part%n, part%n, SPH_D))
        allocate(part%div_nor(part%n), part%nor(part%n, d_Omega), part%centre(d_Omega), part%color(part%n))
        part%centre = centre
        part%rayon = rayon
        part%w = dx**d_Omega
        part%R = SPH_I * dx

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

        open(unit = 10, file = nom_fichier // "_scale.dat")
        write (10, *) centre - rayon
        write (10, *) centre + rayon
        close(10)
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! initialise la pression pour toutes les particules
    ! -------------------------------------------------------------------------------------------------------
    ! part : liste de particules
    subroutine init_var_bulle(centreBulle, rayonBulle, part, func_isInside)
        ! paramètres
        real(rp), dimension(:), intent(in) :: centreBulle
        real(rp), intent(in) :: rayonBulle
        type(Particules), intent(inout) :: part
        logical, external :: func_isInside

        ! variables locales
        integer :: i

        part%rho = 0.0_rp
        part%u = 0.0_rp
        part%P = 2.0_rp * DONNEES_SIGMA / rayonBulle
        call set_dWij(part)
        call set_gradR(part, func_isInside)
    end subroutine

    subroutine init_var_pave(bornes, part, func_isInside)
        ! paramètres
        real(rp), dimension(:, :), intent(in) :: bornes
        type(Particules), intent(inout) :: part
        logical, external :: func_isInside

        ! variables locales
        integer :: i

        part%rho = 0.0_rp
        part%u = 0.0_rp
        part%P = 0.0_rp
        call set_dWij(part)
        !TODO à corriger !!!
        call set_gradR(part, func_isInside)
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! enregistres variables rho u P (éq d'Euler)
    ! -------------------------------------------------------------------------------------------------------
    subroutine write_var(part, f_rho, f_u, f_P)
        ! paramètres
        type(Particules), intent(in) :: part
        character(len=*), intent(in) :: f_rho, f_u, f_P

        ! variables locales
        integer :: i

        open(unit = 10, file = f_rho)
        open(unit = 20, file = f_u)
        open(unit = 30, file = f_P)
        do i = 1, part%n
            write (10, *) part%x(i, :), part%rho(i)
            write (20, *) part%x(i, :), part%u(i, :)
            write (30, *) part%x(i, :), part%P(i)
        end do
        close(10)
        close(20)
        close(30)
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
        real(rp), dimension(SPH_D) :: grad
        real(rp), dimension(part%n, SPH_D * 2) :: nvec
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



    ! -------------------------------------------------------------------------------------------------------
    ! somme FTS sur quartiers
    ! -------------------------------------------------------------------------------------------------------
    subroutine quarter(part, centre)
        ! paramètres
        type(Particules), intent(in) :: part
        real(rp), dimension(2), intent(in) :: centre

        ! variables locales
        real(rp), dimension(4, SPH_D) :: somme, fts, wGP
        integer :: i
        real(rp) :: temp
        real(rp), dimension(SPH_D) :: grad_pressure

        somme = 0.0_rp
        fts = 0.0_rp
        wGP = 0.0_rp
        temp = 0.5_rp * part%dx

        !do i = 1, part%n
        !    if ((part%x(i, 1) >= centre(1) - temp) .and. (part%x(i, 2) >= centre(2) - temp)) then
        !        somme(1, :) = somme(1, :) + part%fts(i, :)
        !    end if
        !    if ((part%x(i, 1) <= centre(1) + temp) .and. (part%x(i, 2) >= centre(2) - temp)) then
        !        somme(2, :) = somme(2, :) + part%fts(i, :)
        !    end if
        !    if ((part%x(i, 1) <= centre(1) + temp) .and. (part%x(i, 2) <= centre(2) + temp)) then
        !        somme(3, :) = somme(3, :) + part%fts(i, :)
        !    end if
        !    if ((part%x(i, 1) >= centre(1) - temp) .and. (part%x(i, 2) <= centre(2) + temp)) then
        !        somme(4, :) = somme(4, :) + part%fts(i, :)
        !    end if
        !end do

        do i = 1, part%n
            if ((part%x(i, 1) >= centre(1) - temp) .and. (part%x(i, 2) >= centre(2) - temp)) then
                call GR_p(i, part, part%P, grad_pressure)
                fts(1, :) = fts(1, :) + part%fts(i, :)
                wGP(1, :) = wGP(1, :) + part%w(i) * grad_pressure
                somme(1, :) = somme(1, :) - part%w(i) * grad_pressure + part%fts(i, :)
            end if
            if ((part%x(i, 1) <= centre(1) + temp) .and. (part%x(i, 2) >= centre(2) - temp)) then
                call GR_p(i, part, part%P, grad_pressure)
                fts(2, :) = fts(2, :) + part%fts(i, :)
                wGP(2, :) = wGP(2, :) + part%w(i) * grad_pressure
                somme(2, :) = somme(2, :) - part%w(i) * grad_pressure + part%fts(i, :)
            end if
            if ((part%x(i, 1) <= centre(1) + temp) .and. (part%x(i, 2) <= centre(2) + temp)) then
                call GR_p(i, part, part%P, grad_pressure)
                fts(3, :) = fts(3, :) + part%fts(i, :)
                wGP(3, :) = wGP(3, :) + part%w(i) * grad_pressure
                somme(3, :) = somme(3, :) - part%w(i) * grad_pressure + part%fts(i, :)
            end if
            if ((part%x(i, 1) >= centre(1) - temp) .and. (part%x(i, 2) <= centre(2) + temp)) then
                call GR_p(i, part, part%P, grad_pressure)
                fts(4, :) = fts(4, :) + part%fts(i, :)
                wGP(4, :) = wGP(4, :) + part%w(i) * grad_pressure
                somme(4, :) = somme(4, :) - part%w(i) * grad_pressure + part%fts(i, :)
            end if
        end do

        print *
        print *, ">>>>>>>>>>>>"
        print *, "Somme des vecteurs D rho w u / Dt par quartier"
        print *, "__________________________________________________________________________________________&
            &_____________________________________"
        print *, "                                                                  |"
        print *, "fts         ", fts(2, :),   " | ", "fts         ", fts(1, :)
        print *, "wGP         ", wGP(2, :),   " | ", "wGP         ", wGP(1, :)
        print *, "top left    ", somme(2, :), " | ", "top right   ", somme(1, :)
        print *, "__________________________________________________________________|_______________________&
            &_____________________________________"
        print *, "                                                                  |"
        print *, "fts         ", fts(3, :),   " | ", "fts         ", fts(4, :)
        print *, "wGP         ", wGP(3, :),   " | ", "wGP         ", wGP(4, :)
        print *, "bottom left ", somme(3, :), " | ", "bottom right", somme(4, :)
        print *, "__________________________________________________________________|_______________________&
            &_____________________________________"
        print *
        print *, "total       ", somme(1, :) + somme(2, :) + somme(3, :) + somme(4, :)
        print *, "<<<<<<<<<<<<"
        print *
    end subroutine

END MODULE donnees
