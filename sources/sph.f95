! ===========================================================================================================
! ===========================================================================================================

MODULE sph

    use math

    implicit none

contains

    ! -------------------------------------------------------------------------------------------------------
    ! créer un quadrillage à partir d'une subdivision de x et de y
    ! -------------------------------------------------------------------------------------------------------
    ! x1 : subdivision des abcisse
    ! x2 : subdivision des ordonnées
    subroutine meshgrid(x1, x2, x)
        ! paramètres
        real(rp), dimension(:), intent(in) :: x1, x2
        real(rp), dimension(:, :), allocatable, intent(out) :: x

        ! variables locales
        integer :: i, j, k

        allocate(x(size(x1) * size(x2), 2))
        k = 1

        do i = 1, size(x1)
            do j = 1, size(x2)
                x(k, :) = (/ x1(i), x2(j) /)
                k = k + 1
            end do
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Maille un disque (plus grand disque centré inclu dans x1 * x2)
    ! -------------------------------------------------------------------------------------------------------
    ! x1 : subdivision abscisse
    ! x2 : subdivision ordonées
    ! x : données du maillage retourné
    ! rayonBulle : rayon du disque
    subroutine meshgridCircle(x1, x2, x, rayonBulle, centreBulle)
        ! paramètres
        real(rp), dimension(:), intent(in) :: x1, x2
        real(rp), dimension(:, :), allocatable, intent(out) :: x
        real(rp), intent(out) :: rayonBulle
        real(rp), dimension(2), intent(out) :: centreBulle

        ! variables locales
        integer :: i, j, k
        real(rp), dimension(size(x1) * size(x2), 2) :: xtemp
        real(rp), dimension(2) :: temp
        real(rp) :: delta

        delta = (x1(2) - x1(1)) / 2.0_rp

        centreBulle = (/ (x1(size(x1)) - x1(1)) / 2.0_rp, (x2(size(x2)) - x2(1)) / 2.0_rp /)
        rayonBulle = maxval(centreBulle) + delta
        centreBulle = centreBulle + (/ x1(1), x2(1) /)

        k = 0

        do i = 1, size(x1)
            do j = 1, size(x2)
                temp = (/ x1(i), x2(j) /)
                if (fnorme2(temp - centreBulle) <= rayonBulle - delta) then
                    k = k + 1
                    xtemp(k, :) = temp
                end if
            end do
        end do

        allocate(x(k, 2))
        x = xtemp(1:k, :)
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Maille un disque étant donné son centre, son rayon et le nombre de particules sur le diamètre
    ! -------------------------------------------------------------------------------------------------------
    subroutine meshCircle(centre, rayon, n_diametre, x)
        ! paramètres
        real(rp), dimension(2), intent(in) :: centre
        real(rp), intent(in) :: rayon
        integer, intent(in) :: n_diametre
        real(rp), dimension(:, :), allocatable, intent(out) :: x

        ! variables locales
        integer :: i, j, k
        real(rp), dimension(n_diametre + 1) :: xm1, xm2
        real(rp), dimension(n_diametre) :: x1, x2
        real(rp), dimension(n_diametre * n_diametre, 2) :: xtemp
        real(rp), dimension(2) :: temp
        real(rp) :: delta

        xm1 = linspace(centre(1) - rayon, centre(1) + rayon, n_diametre + 1)
        xm2 = linspace(centre(2) - rayon, centre(2) + rayon, n_diametre + 1)

        do i = 1, n_diametre
            x1(i) = (xm1(i) + xm1(i + 1)) / 2.0_rp
            x2(i) = (xm2(i) + xm2(i + 1)) / 2.0_rp
        end do

        delta = (x1(2) - x1(1)) / 2.0_rp
        k = 0
        do i = 1, n_diametre
            do j = 1, n_diametre
                temp = (/ x1(i), x2(j) /)
                if (fnorme2(temp - centre) <= rayon - 0.99_rp * delta) then
                    k = k + 1
                    xtemp(k, :) = temp
                end if
            end do
        end do

        print *, k
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

        ! variables locales
        real(rp) :: C

        C = 7.0_rp

        if ((0.0_rp <= q) .and. (q <= 1.0_rp)) then
            theta = C * (1.0_rp - q)**4 * (1.0_rp + 4.0_rp*q)
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
        real(rp), dimension(2), intent(in) :: x
        real(rp), intent(in) :: R

        ! return
        real(rp) :: W

        ! variables locales
        real(rp) :: norm

        call norme2(x, norm)
        if (norm <= R) then
            W = theta(norm / R) / (pi * R**2)
        else
            W = 0.0_rp
        end if
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Gradient en x de W(x, y, R) = W(x - y, R)
    ! -------------------------------------------------------------------------------------------------------
    ! x : élt de R**2 (en pratique coordonnées d'une particule)
    ! R : rayon SPH
    ! grad : résultat (R**2)
    subroutine dx_W_SPH(z, R, grad)
        ! paramètres
        real(rp), dimension(2), intent(in) :: z
        real(rp), intent(in) :: R
        real(rp), dimension(2), intent(out) :: grad

        ! variables locales
        !real(rp) :: u, v
        !real(rp), dimension(2) :: up, vp
        real(rp) :: norm
        real(rp) :: q

        call norme2(z, norm)

        if (norm <= R) then
            q = norm / R
            grad = (7.0_rp / (pi * R**2)) * (-20.0_rp / R) * (z / norm) * q * (1-q)**3
        else
            grad = (/ 0.0_rp, 0.0_rp /)
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Approximation régularisée    Permet vérification Somme(wj * W_SPH_ij) = 1 pour f = 1 partout
    ! -------------------------------------------------------------------------------------------------------
    ! i : numéro de la particule
    ! x : tableau de l'ensemble des particules (une particule de 2 coordonnées par ligne)
    ! w : tableau de l'ensemble des volumes associé aux particules
    ! R : rayon SPH
    ! f : f(x_j) j = 1, np (nombre particules)
    ! image : approx de f(x_i) retournée
    subroutine AR(i, x, w, R, f, image)
        ! paramètres
        integer, intent(in) :: i
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), intent(in) :: R
        real(rp), dimension(:), intent(in) :: f
        real(rp), intent(out) :: image

        ! variables locales
        integer :: j

        image = 0.0_rp

        do j = 1, size(w)
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                image = image + w(j) * f(j) * W_SPH(x(i, :) - x(j, :), R)
            end if
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Gradient régularisé      Permet vérification Somme(wj * \nabla_x W_SPH_ij) = 0 pour f = 1 partout
    ! Trois approximations différentes du gradient de f en x_i (GR, GR_m et GR_p)
    ! -------------------------------------------------------------------------------------------------------
    ! i : numéro de la particule
    ! x : tableau de l'ensemble des particules (une particule de 2 coordonnées par ligne)
    ! w : tableau de l'ensemble des volumes associé aux particules
    ! R : rayon SPH
    ! f : f(x_j) j = 1, np (nombre particules)
    ! image : approx de \nabla f(x_i) retournée
    subroutine GR(i, x, w, R, f, image)
        ! paramètres
        integer, intent(in) :: i
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), intent(in) :: R
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(2), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(2) :: grad

        image = 0.0_rp

        do j = 1, i - 1
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * f(j) * grad
            end if
        end do

        do j = i + 1, size(w)
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * f(j) * grad
            end if
        end do
    end subroutine

    ! i : numéro de la particule
    ! x : tableau de l'ensemble des particules (une particule de 2 coordonnées par ligne)
    ! w : tableau de l'ensemble des volumes associé aux particules
    ! R : rayon SPH
    ! f : f(x_j) j = 1, np (nombre particules)
    ! image : approx de \nabla f(x_i) retournée
    subroutine GR_m(i, x, w, R, f, image)
        ! paramètres
        integer, intent(in) :: i
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), intent(in) :: R
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(2), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(2) :: grad

        image = 0.0_rp

        do j = 1, i - 1
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * (f(j) - f(i)) * grad
            end if
        end do

        do j = i + 1, size(w)
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * (f(j) - f(i)) * grad
            end if
        end do
    end subroutine

    ! i : numéro de la particule
    ! x : tableau de l'ensemble des particules (une particule de 2 coordonnées par ligne)
    ! w : tableau de l'ensemble des volumes associé aux particules
    ! R : rayon SPH
    ! f : f(x_j) j = 1, np (nombre particules)
    ! image : approx de \nabla f(x_i) retournée
    subroutine GR_p(i, x, w, R, f, image)
        ! paramètres
        integer, intent(in) :: i
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), intent(in) :: R
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(2), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(2) :: grad

        image = 0.0_rp

        do j = 1, i - 1
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * (f(j) + f(i)) * grad
            end if
        end do

        do j = i + 1, size(w)
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * (f(j) + f(i)) * grad
            end if
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Noyau cohésion pour tension de surface (cf C(r) Akinci p. 3)
    ! -------------------------------------------------------------------------------------------------------
    function C_Akinci(r, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: r, R_SPH

        ! return
        real(rp) :: C

        ! variables locales
        real(rp) :: kernel_TS, q
        real(RP), save :: TSkernorm = 139.0_rp / 1120.0_rp * 336.0_rp / 37.0_rp
        ! en dimension 2, sinon vaut 1 en dimension 3
        real(RP), save :: C_norm = 417.0_rp / 370.0_rp
        integer :: d

        ! dimension
        !D = 2
        !q = r/R_SPH

        !kernel_TS = 0._rp
        !if ((2._rp*q>1._rp).and.(q<=1._RP)) then
        !    kernel_TS = (1._rp-q)**3._rp*q**3._rp
        !elseif ((q>0._rp).and.(2._rp*q<=1._rp)) then
        !    kernel_TS = 2._rp*(1._rp-q)**3._rp*q**3._rp - 1._rp/64._rp
        !endif
        !kernel_TS = TSkernorm*32._RP/(PI*R_SPH**real(D, rp))*kernel_TS

        !if ((0.0_rp < r) .and. (r <= R_SPH / 2.0_rp)) then
        !    C = (32.0_rp / (pi * R_SPH**9)) * (2.0_rp * (R_SPH - r)**3 * r**3 - R_SPH**6 / 64.0_rp)
        !else if ((R_SPH / 2.0_rp < r) .and. (r <= R_SPH)) then
        !    C = (32.0_rp / (pi * R_SPH**9)) * (R_SPH - r)**3 * r**3
        !else
        !    C = 0.0_rp
        !end if

        d = 2
        if ((0.0_rp < r) .and. (r <= R_SPH / 2.0_rp)) then
            C = C_norm * (32.0_rp / (pi * R_SPH**d * R_SPH**6)) * (2.0_rp * (R_SPH - r)**3 * r**3 - R_SPH**6 / 64.0_rp)
        else if ((R_SPH / 2.0_rp < r) .and. (r <= R_SPH)) then
            C = C_norm * (32.0_rp / (pi * R_SPH**d * R_SPH**6)) * (R_SPH - r)**3 * r**3
        else
            C = 0.0_rp
        end if
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! vecteur ni = R_SPH * \sum_j ( wj \nabla W_ij ) pour toute particule i
    ! -------------------------------------------------------------------------------------------------------
    ! R_SPH : rayon du noyau SPH
    ! x : particules
    ! w : volume des particules
    ! n : vecteurs normaux non normalisés (en sortie)
    subroutine normale(R_SPH, x, w, n)
        ! paramètres
        real(rp), intent(in) :: R_SPH
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), dimension(size(w), 2), intent(out) :: n

        ! variables locales
        integer :: i, j, k
        real(rp), dimension(2) :: ni, grad

        do i = 1, size(w)
            ni = 0.0_rp
            do j = 1, i - 1
                if (fnorme2(x(i, :) - x(j, :)) <= R_SPH) then
                    call dx_W_SPH(x(i, :) - x(j, :), R_SPH, grad)
                    ni = ni + w(j) * grad
                end if
            end do

            do j = i + 1, size(w)
                if (fnorme2(x(i, :) - x(j, :)) <= R_SPH) then
                    call dx_W_SPH(x(i, :) - x(j, :), R_SPH, grad)
                    ni = ni + w(j) * grad
                end if
            end do

            n(i, :) = R_SPH * ni
        end do

    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Force de tension de surface solver SPH
    ! -------------------------------------------------------------------------------------------------------
    ! y : coefficient de tension de surface (gamma)
    ! i : numéro d'une particule
    ! x : particules (coordonnées)
    ! w : volume des particules
    ! R_SPH : rayon SPH
    subroutine F_TS(y, i, x, w, n, R_SPH, F)
        ! paramètres
        real(rp), intent(in) :: y
        integer, intent(in) :: i
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), dimension(:, :), intent(in) :: n
        real(rp), intent(in) :: R_SPH
        real(rp), dimension(2), intent(out) :: F

        ! variables locales
        real(rp), dimension(2) :: ni, nj
        integer :: k, np, j

        np = size(w)

        F = 0.0_rp
        do j = 1, i - 1
            if (fnorme2(x(i, :) - x(j, :)) <= R_SPH) then
                F = F + y * w(i) * R_SPH * C_Akinci(fnorme2(x(i, :) - x(j, :)), R_SPH) &
                    * (x(i, :) - x(j, :)) / fnorme2(x(i, :) - x(j, :)) &
                    + y * R_SPH * (n(i, :) - n(j, :))
            end if
        end do
        do j = i + 1, np
            if (fnorme2(x(i, :) - x(j, :)) <= R_SPH) then
                F = F + y * w(i) * R_SPH * C_Akinci(fnorme2(x(i, :) - x(j, :)), R_SPH) &
                    * (x(i, :) - x(j, :)) / fnorme2(x(i, :) - x(j, :)) &
                    + y * R_SPH * (n(i, :) - n(j, :))
            end if
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! noyau Kordilla
    ! -------------------------------------------------------------------------------------------------------
    function C_Kordilla(x, R) result(W)
        ! paramètres
        real(rp), dimension(2), intent(in) :: x
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
        real(rp), dimension(2), intent(out) :: F

        ! variables locales
        real(rp) :: s_ff, A, B, h1, h2
        real(rp), dimension(2) :: r
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
    ! x : particules (coordonnées)
    ! w : volume des particules
    ! R_SPH : rayon SPH
    subroutine F_TS_cohesion(y, i, x, w, R_SPH, F)
        ! paramètres
        real(rp), intent(in) :: y
        integer, intent(in) :: i
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), intent(in) :: R_SPH
        real(rp), dimension(2), intent(out) :: F

        ! variables locales
        integer :: k, np, j

        np = size(w)

        F = 0.0_rp
        do j = 1, i - 1
            if (fnorme2(x(i, :) - x(j, :)) <= R_SPH) then
                F = F + y * w(i) * R_SPH * C_Akinci(fnorme2(x(i, :) - x(j, :)), R_SPH) &
                    * (x(i, :) - x(j, :)) / fnorme2(x(i, :) - x(j, :))
            end if
        end do
        do j = i + 1, np
            if (fnorme2(x(i, :) - x(j, :)) <= R_SPH) then
                F = F + y * w(i) * R_SPH * C_Akinci(fnorme2(x(i, :) - x(j, :)), R_SPH) &
                    * (x(i, :) - x(j, :)) / fnorme2(x(i, :) - x(j, :))
            end if
        end do
    end subroutine

END MODULE sph
