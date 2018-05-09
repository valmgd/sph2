! ===========================================================================================================
! Code pour évaluation de la tension de surface et tests sur schéma SPH sur une initialisation
! (pas d'itérations en temps)
!
! Prévu pour fonctionner en 2D ou 3D
! Pour cela changer la constante d. Les autres variables avec exposant sont des constantes qui dépendent
! de la dimension. Elles sont ainsi calculées une seule fois. L'astuce est la suivante :
!
! a**(3-d) * b**(d-2) = {a si d == 2, b si d == 3}
! ===========================================================================================================

MODULE sph

    use math
    use var

    implicit none

contains

    ! -------------------------------------------------------------------------------------------------------
    ! destructeur
    ! -------------------------------------------------------------------------------------------------------
    ! part : variable structuré de type Particules que l'on veut désallouer
    subroutine rm_Particules(part)
        ! paramètres
        type(Particules), intent(inout) :: part

        if (allocated(part%x)) then
            deallocate(part%x)
        end if
        if (allocated(part%w)) then
            deallocate(part%w)
        end if

        if (allocated(part%gradR)) then
            deallocate(part%gradR)
        end if
        if (allocated(part%fts)) then
            deallocate(part%fts)
        end if

        if (allocated(part%rho)) then
            deallocate(part%rho)
        end if
        if (allocated(part%u)) then
            deallocate(part%u)
        end if
        if (allocated(part%P)) then
            deallocate(part%P)
        end if
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



    ! -------------------------------------------------------------------------------------------------------
    ! fonction pour noyau SPH
    ! -------------------------------------------------------------------------------------------------------
    ! q : q en dehors de [0, 1] <=> x en dehors du support de W
    function theta(q)
        ! paramètres
        real(rp), intent(in) :: q

        ! return
        real(rp) :: theta

        if ((0.0_rp <= q) .and. (q <= 1.0_rp)) then
            theta = SPH_C_NORM_NOYAU * (1.0_rp - q)**4 * (1.0_rp + 4.0_rp*q)
        else
            theta = 0.0_rp
        end if
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Noyau SPH
    ! -------------------------------------------------------------------------------------------------------
    ! x : élt de R**2 (en pratique coordonnées d'une particule)
    ! R : rayon SPH
    ! W : résultat (réel)
    function W_SPH(x, R) result(W)
        ! paramètres
        real(rp), dimension(SPH_D), intent(in) :: x
        real(rp), intent(in) :: R

        ! return
        real(rp) :: W

        ! variables locales
        real(rp) :: norm

        call norme2(x, norm)
        if (norm <= R) then
            W = theta(norm / R) / (SPH_VOLUME_SUPP * R**SPH_D)
        else
            W = 0.0_rp
        end if
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Gradient en x de W(x, y, R) = W(x - y, R)
    ! -------------------------------------------------------------------------------------------------------
    ! z : élt de R**2 (en pratique coordonnées d'une particule)
    ! R : rayon SPH
    ! grad : résultat (R**2)
    subroutine dx_W_SPH(z, R, grad)
        ! paramètres
        real(rp), dimension(SPH_D), intent(in) :: z
        real(rp), intent(in) :: R
        real(rp), dimension(SPH_D), intent(out) :: grad

        ! variables locales
        !real(rp) :: u, v
        !real(rp), dimension(2) :: up, vp
        real(rp) :: norm
        real(rp) :: q

        call norme2(z, norm)

        if (norm <= R) then
            q = norm / R
            grad = (SPH_C_NORM_NOYAU / (SPH_VOLUME_SUPP * R**SPH_D)) * (-20.0_rp / R) * (z / norm) * q * (1-q)**3
        else
            grad = 0.0_rp
        end if
    end subroutine



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
                call dx_W_SPH(part%x(i, :) - part%x(j, :), part%R, grad)
                image = image + part%w(j) * f(j) * grad
            end if
        end do

        do j = i + 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                call dx_W_SPH(part%x(i, :) - part%x(j, :), part%R, grad)
                image = image + part%w(j) * f(j) * grad
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
                call dx_W_SPH(part%x(i, :) - part%x(j, :), part%R, grad)
                image = image + part%w(j) * (f(j) - f(i)) * grad
            end if
        end do

        do j = i + 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                call dx_W_SPH(part%x(i, :) - part%x(j, :), part%R, grad)
                image = image + part%w(j) * (f(j) - f(i)) * grad
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
                call dx_W_SPH(part%x(i, :) - part%x(j, :), part%R, grad)
                image = image + part%w(j) * (f(j) + f(i)) * grad
            end if
        end do

        do j = i + 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                call dx_W_SPH(part%x(i, :) - part%x(j, :), part%R, grad)
                image = image + part%w(j) * (f(j) + f(i)) * grad
            end if
        end do
    end subroutine



    ! =======================================================================================================
    ! TENSION SUPERFICIELLE
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! Noyau cohésion pour tension de surface (cf C(r) Akinci p. 3)
    ! -------------------------------------------------------------------------------------------------------
    function C_Akinci(r, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: r, R_SPH

        ! return
        real(rp) :: C

        if ((0.0_rp < r) .and. (r <= R_SPH / 2.0_rp)) then
            C = SPH_C_NORM_AKINCI * (32.0_rp / (pi * R_SPH**SPH_D * R_SPH**6)) * &
                (2.0_rp * (R_SPH - r)**3 * r**3 - R_SPH**6 / 64.0_rp)

        else if ((R_SPH / 2.0_rp < r) .and. (r <= R_SPH)) then
            C = SPH_C_NORM_AKINCI * (32.0_rp / (pi * R_SPH**SPH_D * R_SPH**6)) * &
                (R_SPH - r)**3 * r**3

        else
            C = 0.0_rp
        end if
    end function



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
                    call dx_W_SPH(part%x(i, :) - part%x(j, :), part%R, grad)
                    ni = ni + part%w(j) * grad
                end if
            end do

            do j = i + 1, part%n
                if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                    call dx_W_SPH(part%x(i, :) - part%x(j, :), part%R, grad)
                    ni = ni + part%w(j) * grad
                end if
            end do

            part%gradR(i, :) = part%R * ni
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Force de tension de surface solver SPH
    ! -------------------------------------------------------------------------------------------------------
    ! sigma : coefficient de tension de surface (gamma)
    ! i : numéro d'une particule
    ! part : liste des particules
    subroutine FTS_akinci(sigma, i, part, F)
        ! paramètres
        real(rp), intent(in) :: sigma
        integer, intent(in) :: i
        type(Particules), intent(in) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        ! variables locales
        real(rp), dimension(SPH_D) :: ni, nj
        integer :: k, j

        F = 0.0_rp
        do j = 1, i - 1
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                F = F + sigma * part%w(i) * part%R * C_Akinci(fnorme2(part%x(i, :) - part%x(j, :)), part%R) &
                    * (part%x(i, :) - part%x(j, :)) / fnorme2(part%x(i, :) - part%x(j, :)) &
                    + sigma * part%R * (part%gradR(i, :) - part%gradR(j, :))
            end if
        end do
        do j = i + 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                F = F + sigma * part%w(i) * part%R * C_Akinci(fnorme2(part%x(i, :) - part%x(j, :)), part%R) &
                    * (part%x(i, :) - part%x(j, :)) / fnorme2(part%x(i, :) - part%x(j, :)) &
                    + sigma * part%R * (part%gradR(i, :) - part%gradR(j, :))
            end if
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! calculer la tension de surface pour toute particule i
    ! -------------------------------------------------------------------------------------------------------
    ! sigma : coefficient de tension superficielle (cste)
    ! part : liste des particules
    subroutine set_fts(FTS_func, sigma, part)
        ! paramètres
        interface
            subroutine FTS_func(sigma_, i_, part_, F)
                use math
                use var
                real(rp), intent(in) :: sigma_
                integer, intent(in) :: i_
                type(Particules), intent(in) :: part_
                real(rp) :: part
                real(rp), dimension(SPH_D), intent(out) :: F
            end subroutine
        end interface
        real(rp), intent(in) :: sigma
        type(Particules), intent(inout) :: part

        ! variables locales
        integer :: i

        do i = 1, part%n
            call FTS_func(sigma, i, part, part%fts(i, :))
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! noyau Kordilla
    ! -------------------------------------------------------------------------------------------------------
    function C_Kordilla(x, R) result(W)
        ! paramètres
        real(rp), dimension(SPH_D), intent(in) :: x
        real(rp), intent(in) :: R

        ! return
        real(rp) :: W

        ! variables locales
        real(rp) :: length

        length = fnorme2(x)

        if ((0.0_rp <= length) .and. (length < R / 3.0_rp)) then
            W = (3.0_rp - 3.0_rp * length / R)**5 - 6.0_rp * (2.0_rp - 3.0_rp * length / R)**5 &
                + 15.0_rp * (1.0_rp - 3.0_rp * length / R)**5
            W = W * (81.0_rp / (359.0_rp * pi * R**3.0_rp))
        else if ((R / 3.0_rp < length) .and. (length < 2.0_rp * R / 3.0_rp)) then
            W = (3.0_rp - 3.0_rp * length / R)**5 - 6.0_rp * (2.0_rp - 3.0_rp * length / R)**5
            W = W * (81.0_rp / (359.0_rp * pi * R**3.0_rp))
        else if ((2.0_rp * R / 3.0_rp <= length) .and. (length < R)) then
            W = (3.0_rp - 3.0_rp * length / R)**5
            W = W * (81.0_rp / (359.0_rp * pi * R**3.0_rp))
        else
            W = 0.0_rp
        end if
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Force d'une particule sur une autre (Kordilla)
    ! -------------------------------------------------------------------------------------------------------
    ! y : coeff de tension de surface (gamma)
    subroutine F_pp_kordilla(y, R_SPH, i, j, x, F)
        ! paramètres
        real(rp), intent(in) :: y, R_SPH
        integer, intent(in) :: i, j
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(SPH_D), intent(out) :: F

        ! variables locales
        real(rp) :: s_ff, A, B, h1, h2
        real(rp), dimension(SPH_D) :: r
        real(rp) :: length

        A = 2.0_rp
        B = -1.0_rp
        h1 = 0.8_rp
        h2 = 1.0_rp

        r = x(i, :) - x(j, :)
        length = fnorme2(r)

        if (length <= R_SPH) then
            F = y * (A * C_Kordilla(x(i, :), h1) * r / length + B * C_Kordilla(x(i, :), h2) * r / length)
        else
            F = 0.0_rp
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! force de tension de surface Kordilla
    ! -------------------------------------------------------------------------------------------------------
    subroutine F_TS_kordilla()
        ! paramètres
        !++!

        ! variables locales
        !++!
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Force de tension de surface idem akinci mais avec juste F_cohesion
    ! -------------------------------------------------------------------------------------------------------
    ! y : coefficient de tension de surface (gamma)
    ! i : numéro d'une particule
    ! part : liste des particules
    subroutine F_TS_cohesion(y, i, part, F)
        ! paramètres
        real(rp), intent(in) :: y
        integer, intent(in) :: i
        type(Particules), intent(in) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        ! variables locales
        integer :: k, j

        F = 0.0_rp
        do j = 1, i - 1
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                F = F + y * part%w(i) * part%R * C_Akinci(fnorme2(part%x(i, :) - part%x(j, :)), part%R) &
                    * (part%x(i, :) - part%x(j, :)) / fnorme2(part%x(i, :) - part%x(j, :))
            end if
        end do
        do j = i + 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                F = F + y * part%w(i) * part%R * C_Akinci(fnorme2(part%x(i, :) - part%x(j, :)), part%R) &
                    * (part%x(i, :) - part%x(j, :)) / fnorme2(part%x(i, :) - part%x(j, :))
            end if
        end do
    end subroutine



    ! =======================================================================================================
    ! Schéma SPH
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! ressemble à une itération en temps de l'une des équation du schéma SPH
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

            call prodScal(part%x(i, :) - centre, d_rwu_dt(i, :), prod)
            plot_vec(i, :) = (/ part%x(i, :), d_rwu_dt(i, :), prod /)
        end do
        !print *, "sum_i [ w_i GR_p(P)_i + (FTS)_i ] =", sum(d_rwu_dt(:, 1)), sum(d_rwu_dt(:, 2))
        !print *, sum(d_rwu_dt(:, 1)), sum(d_rwu_dt(:, 2))
        !call affMat(d_rwu_dt)

        call writeMat(plot_vec, "../sorties/plot_vec.dat")
    end subroutine

END MODULE sph
