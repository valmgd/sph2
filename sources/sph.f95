! ===========================================================================================================
! Code pour évaluation de la tension de surface et tests sur schéma SPH sur une initialisation
! (pas d'itérations en temps)
!
! Prévu pour fonctionner en 2D ou 3D
! Pour cela changer la constante d. Les autres variables avec exposant sont des constantes qui dépendent
! de la dimension. Elles sont ainsi calculées une seule fois. L'astuce est la suivante :
!
! a**(3-d) * b**(d-2) = {a si d == 2, b si d == 3}
!
!   subroutine set_gradR(part)
!   subroutine meshgrid(xyz, x)
!   subroutine meshCircle(d_Omega, centre, rayon, n_diametre, x)
!   subroutine dx_W_SPH(z, R, grad)
!   subroutine AR(i, part, f, image)
!   subroutine GR(i, part, f, image)
!   subroutine GR_m(i, part, f, image)
!   subroutine GR_p(i, part, f, image)
!   subroutine DR(i, part, f, image)
!   subroutine DR_m(i, part, f, image)
!   subroutine DR_p(i, part, f, image)
!   subroutine iter_SPH(part, centre)
! ===========================================================================================================

MODULE sph

    use tsup

    implicit none

contains

    ! =======================================================================================================
    ! MAILLAGE
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! gradient : \nabla W_SPH(xi - xj)
    ! -------------------------------------------------------------------------------------------------------
    ! part : liste des particules
    subroutine set_dWij(part)
        ! paramètres
        type(Particules), intent(inout) :: part

        ! variables locales
        integer :: i, j

        part%dWij = 0.0_rp

        do j = 1, part%n
            do i = 1, j - 1
                if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                    call dx_W_SPH(part%x(i, :) - part%x(j, :), part%R, part%dWij(i, j, :))
                end if
            end do
            do i = j + 1, part%n
                if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                    call dx_W_SPH(part%x(i, :) - part%x(j, :), part%R, part%dWij(i, j, :))
                end if
            end do
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! vecteur ni = R_SPH * \sum_j ( wj \nabla W_ij ) pour toute particule i
    ! -------------------------------------------------------------------------------------------------------
    ! part : liste des particules
    ! n : vecteurs normaux non normalisés (en sortie)
    subroutine set_gradR(part)
        ! paramètres
        type(Particules), intent(inout) :: part

        ! variables locales
        integer :: i, j, k
        real(rp), dimension(SPH_D) :: ni, grad

        do i = 1, part%n
            ni = 0.0_rp
            do j = 1, i - 1
                if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                    ni = ni + part%w(j) * part%dWij(i, j, :)
                end if
            end do

            do j = i + 1, part%n
                if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                    ni = ni + part%w(j) * part%dWij(i, j, :)
                end if
            end do

            part%gradR(i, :) = part%R * ni
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! créer un quadrillage à partir d'une subdivision de x1, de x2 et éventuellement de x3
    ! -------------------------------------------------------------------------------------------------------
    ! xyz : 2 ou 3 subdivisions de chacun des 2 ou 3 axes d'espace (en colonne)
    ! x : maillage 2 ou 3D retourné
    subroutine meshgrid(xyz, x)
        ! paramètres
        real(rp), dimension(:, :), intent(in) :: xyz
        real(rp), dimension(:, :), allocatable, intent(out) :: x

        ! variables locales
        integer :: i, j
        integer :: d_dom, length, N
        integer, dimension(:), allocatable :: indice

        length = size(xyz, 1)
        d_dom = size(xyz, 2)
        N = length**d_dom

        allocate(indice(d_dom), x(N, d_dom))
        indice = 1

        do i = 1, N
            do j = 1, d_dom
                x(i, j) = xyz(indice(j), j)
            end do

            do j = d_dom, 1, -1
                if (indice(j) < length) then
                    indice(j) = indice(j) + 1
                    exit
                else
                    indice(j) = 1
                end if
            end do
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Maille un disque étant donné son centre, son rayon et le nombre de particules sur le diamètre
    ! -------------------------------------------------------------------------------------------------------
    ! centre : centre de la bulle à mailler
    ! rayon : rayon de la bulle à mailler
    ! n_diametre : nombre de points voulus sur un diametre
    ! x : coordonnées des particules retournées
    subroutine meshCircle(d_Omega, centre, rayon, n_diametre, x)
        ! paramètres
        integer, intent(in) :: d_Omega
        real(rp), dimension(d_Omega), intent(in) :: centre
        real(rp), intent(in) :: rayon
        integer, intent(in) :: n_diametre
        real(rp), dimension(:, :), allocatable, intent(out) :: x

        ! variables locales
        integer :: i, j, k
        real(rp), dimension(n_diametre + 1, d_Omega) :: xm
        real(rp), dimension(n_diametre, d_Omega) :: axes
        real(rp), dimension(n_diametre**d_Omega, d_Omega) :: xtemp
        real(rp), dimension(d_Omega) :: temp
        real(rp) :: delta
        integer :: length, N
        integer, dimension(:), allocatable :: indice


        do i = 1, d_Omega
            xm(:, i) = linspace(centre(i) - rayon, centre(i) + rayon, n_diametre + 1)
        end do

        do i = 1, n_diametre
            axes(i, :) = (xm(i, :) + xm(i + 1, :)) / 2.0_rp
        end do

        delta = (axes(2, 1) - axes(1, 1)) / 2.0_rp
        length = size(axes, 1)
        N = length**d_Omega

        allocate(indice(d_Omega))
        indice = 1

        k = 0
        do i = 1, N
            do j = 1, d_Omega
                temp(j) = axes(indice(j), j)
            end do
            if (fnorme2(temp - centre) <= rayon - 0.99_rp * delta) then
                k = k + 1
                xtemp(k, :) = temp
            end if

            do j = d_Omega, 1, -1
                if (indice(j) < length) then
                    indice(j) = indice(j) + 1
                    exit
                else
                    indice(j) = 1
                end if
            end do
        end do

        allocate(x(k, 2))
        x = xtemp(1:k, :)
    end subroutine



    ! =======================================================================================================
    ! OPÉRATEURS RÉGULARISÉS
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! Approximation régularisée    Permet vérification Somme(wj * W_SPH_ij) = 1 pour f = 1 partout
    ! -------------------------------------------------------------------------------------------------------
    ! i : numéro de la particule
    ! part : liste des particules
    ! f : f(x_j) j = 1, np (nombre particules)
    ! image : approx de f(x_i) retournée
    subroutine AR(i, part, f, image)
        ! paramètres
        integer, intent(in) :: i
        type(Particules), intent(in) :: part
        real(rp), dimension(:), intent(in) :: f
        real(rp), intent(out) :: image

        ! variables locales
        integer :: j

        image = 0.0_rp

        do j = 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                image = image + part%w(j) * f(j) * W_SPH(part%x(i, :) - part%x(j, :), part%R)
            end if
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Gradient régularisé      Permet vérification Somme(wj * \nabla_x W_SPH_ij) = 0 pour f = 1 partout
    ! Trois approximations différentes du gradient de f en x_i (GR, GR_m et GR_p)
    ! -------------------------------------------------------------------------------------------------------
    ! i : numéro de la particule
    ! part : liste des particules
    ! f : f(x_j) j = 1, np (nombre particules)
    ! image : approx de \nabla f(x_i) retournée
    subroutine GR(i, part, f, image)
        ! paramètres
        integer, intent(in) :: i
        type(Particules), intent(in) :: part
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(SPH_D), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(SPH_D) :: grad

        image = 0.0_rp

        do j = 1, i - 1
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                image = image + part%w(j) * f(j) * part%dWij(i, j, :)
            end if
        end do

        do j = i + 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                image = image + part%w(j) * f(j) * part%dWij(i, j, :)
            end if
        end do
    end subroutine

    ! i : numéro de la particule
    ! part : liste des particules
    ! f : f(x_j) j = 1, np (nombre particules)
    ! image : approx de \nabla f(x_i) retournée
    subroutine GR_m(i, part, f, image)
        ! paramètres
        integer, intent(in) :: i
        type(Particules), intent(in) :: part
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(SPH_D), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(SPH_D) :: grad

        image = 0.0_rp

        do j = 1, i - 1
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                image = image + part%w(j) * (f(j) - f(i)) * part%dWij(i, j, :)
            end if
        end do

        do j = i + 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                image = image + part%w(j) * (f(j) - f(i)) * part%dWij(i, j, :)
            end if
        end do
    end subroutine

    ! i : numéro de la particule
    ! part : liste des particules
    ! f : f(x_j) j = 1, np (nombre particules)
    ! image : approx de \nabla f(x_i) retournée
    subroutine GR_p(i, part, f, image)
        ! paramètres
        integer, intent(in) :: i
        type(Particules), intent(in) :: part
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(SPH_D), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(SPH_D) :: grad

        image = 0.0_rp

        do j = 1, i - 1
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                image = image + part%w(j) * (f(j) + f(i)) * part%dWij(i, j, :)
            end if
        end do

        do j = i + 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                image = image + part%w(j) * (f(j) + f(i)) * part%dWij(i, j, :)
            end if
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Divergence régularisée
    ! Trois approximations différentes de la divergence de f en x_i (DR, DR_m et DR_p)
    ! -------------------------------------------------------------------------------------------------------
    subroutine DR(i, part, f, image)
        ! paramètres
        integer, intent(in) :: i
        type(Particules), intent(in) :: part
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(SPH_D), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(SPH_D) :: grad
    end subroutine

    subroutine DR_m(i, part, f, image)
        ! paramètres
        integer, intent(in) :: i
        type(Particules), intent(in) :: part
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(SPH_D), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(SPH_D) :: grad
    end subroutine

    subroutine DR_p(i, part, f, image)
        ! paramètres
        integer, intent(in) :: i
        type(Particules), intent(in) :: part
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(SPH_D), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(SPH_D) :: grad
    end subroutine



    ! =======================================================================================================
    ! SCHÉMA SPH
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! une itération du schéma SPH
    ! -------------------------------------------------------------------------------------------------------
    subroutine iter_SPH(part, centre)
        ! paramètres
        type(Particules), intent(in) :: part
        real(rp), dimension(SPH_D), intent(in) :: centre

        ! variables locales
        real(rp), dimension(SPH_D) :: grad_pressure
        real(rp), dimension(part%n, SPH_D) :: d_rwu_dt
        integer :: i
        real(rp), dimension(part%n, 2 * SPH_D + 1) :: plot_vec
        real(rp) :: prod

        do i = 1, part%n
            call GR_p(i, part, part%P, grad_pressure)
            d_rwu_dt(i, :) = part%fts(i, :) - part%w(i) * grad_pressure
            ! print *, part%fts(i, :), part%w(i) * grad_pressure

            !call prodScal(part%x(i, :) - centre, d_rwu_dt(i, :), prod)
            !plot_vec(i, :) = (/ part%x(i, :), d_rwu_dt(i, :), prod /)

            !call prodScal(part%x(i, :) - centre, part%w(i) * grad_pressure, prod)
            !plot_vec(i, :) = (/ part%x(i, :), part%w(i) * grad_pressure, prod /)

            call prodScal(part%x(i, :) - centre, part%fts(i, :), prod)
            plot_vec(i, :) = (/ part%x(i, :), part%fts(i, :), prod /)
        end do
        !print *, "sum_i [ -w_i GR_p(P)_i + (FTS)_i ] =", sum(d_rwu_dt(:, 1)), sum(d_rwu_dt(:, 2))
        !print *, sum(d_rwu_dt(:, 1)), sum(d_rwu_dt(:, 2))
        !call affMat(d_rwu_dt)

        call writeMat(plot_vec, "../sorties/plot_vec.dat")
    end subroutine

END MODULE sph
