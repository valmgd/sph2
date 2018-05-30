! ===========================================================================================================
! Module pour tension superficielle
!
!   function C_akinci(r, R_SPH) result(C)
!   subroutine FTS_akinci(sigma, i, part, F)
!   subroutine F_pp_kordilla(y, R_SPH, i, j, x, F)
!   subroutine F_TS_kordilla()
!   subroutine F_TS_cohesion(y, i, part, F)
!   subroutine set_fts(FTS_func, sigma, part)
! ===========================================================================================================

MODULE tsup

    use kernel

    implicit none

contains

    ! =======================================================================================================
    ! TENSION SUPERFICIELLE AKINCI
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! Noyau cohésion pour tension de surface (cf C(r) Akinci p. 3)
    ! -------------------------------------------------------------------------------------------------------
    function C_akinci(q, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: q, R_SPH

        ! return
        real(rp) :: C

        if ((0.0_rp < q) .and. (q <= 0.5_rp)) then
            C = (SPH_NORM_COHESION_AKINCI / R_SPH**SPH_D) * &
                (2.0_rp * (1.0_rp - q)**3 * q**3 - 1.0_rp / 64.0_rp)

        else if ((0.5_rp < q) .and. (q <= 1.0_rp)) then
            C = (SPH_NORM_COHESION_AKINCI / R_SPH**SPH_D) * (1.0_rp - q)**3 * q**3

        else
            C = 0.0_rp
        end if
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Force de cohésion des particules i <- j
    ! -------------------------------------------------------------------------------------------------------
    subroutine F_ij_akinci(sigma, i, j, part, F)
        ! paramètres
        real(rp), intent(in) :: sigma
        integer, intent(in) :: i, j
        type(Particules), intent(in) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
            ! cohesion force
            !TODO : vérifier si il y a un moins devant le sigma
            F = -sigma * part%w(i) * part%R * C_akinci(fnorme2(part%x(i, :) - part%x(j, :)) / part%R, part%R) &
                * (part%x(i, :) - part%x(j, :)) / fnorme2(part%x(i, :) - part%x(j, :))

            ! curvature force
            F = F - sigma * part%R * (part%gradR(i, :) - part%gradR(j, :))
        else
            F = 0.0_rp
        end if
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
        type(Particules), intent(inout) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        ! variables locales
        real(rp), dimension(SPH_D) :: ni, nj, Fij
        integer :: k, j

        F = 0.0_rp
        do j = 1, i - 1
            call F_ij_akinci(sigma, i, j, part, Fij)
            F = F + Fij
        end do
        do j = i + 1, part%n
            call F_ij_akinci(sigma, i, j, part, Fij)
            F = F + Fij
        end do
    end subroutine



    ! =======================================================================================================
    ! TENSION SUPERFICIELLE LIU
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! noyau de cohésion Liu
    ! -------------------------------------------------------------------------------------------------------
    function C_liu(q, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: q, R_SPH

        ! return
        real(rp) :: C

        ! variables locales
        real(rp), dimension(SPH_D) :: F
        real(rp) :: s_ff, A, B, h1, h2
        real(rp), dimension(SPH_D) :: r
        real(rp) :: length

        A = 2.0_rp
        B = 1.0_rp
        h1 = 0.8_rp * R_SPH
        h2 = 1.0_rp * R_SPH

        r = (/ q * R_SPH, 0.0_rp /)
        length = fnorme2(r)

        if (q <= 1.0_rp) then
            C = A * W_SPH_liu(r, h1) - B * W_SPH_liu(r, h2)
        else
            C = 0.0_rp
        end if

    end function

    ! -------------------------------------------------------------------------------------------------------
    ! Force de cohésion d'une particule sur une autre (Liu, puis repris par Kordilla)
    ! -------------------------------------------------------------------------------------------------------
    ! sigma : coeff de tension de surface (gamma)
    subroutine F_ij_liu(sigma, i, j, part, F)
        ! paramètres
        real(rp), intent(in) :: sigma
        integer, intent(in) :: i, j
        type(Particules), intent(in) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        ! variables locales
        real(rp) :: s_ff, A, B, h1, h2
        real(rp), dimension(SPH_D) :: r
        real(rp) :: length

        A = 2.0_rp
        B = -1.0_rp
        h1 = 0.8_rp * part%R
        h2 = 1.0_rp * part%R

        r = part%x(i, :) - part%x(j, :)
        length = fnorme2(r)

        if (length <= part%R) then
            F = sigma * (r / length) * C_liu(length / part%R, part%R)
        else
            F = 0.0_rp
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! force de tension de surface Kordilla
    ! -------------------------------------------------------------------------------------------------------
    subroutine FTS_liu(sigma, i, part, F)
        ! paramètres
        ! paramètres
        real(rp), intent(in) :: sigma
        integer, intent(in) :: i
        type(Particules), intent(inout) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        ! variables locales
        real(rp), dimension(SPH_D) :: Fij
        integer :: k, j

        F = 0.0_rp

        do j = 1, i - 1
            call F_ij_liu(sigma, i, j, part, F)
            F = F + part%w(j) * Fij
        end do
        do j = i + 1, part%n
            call F_ij_liu(sigma, i, j, part, F)
            F = F + part%w(j) * Fij
        end do
    end subroutine



    ! =======================================================================================================
    ! FTS NEW
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! Noyau cohésion pour tension de surface (cf C(r) Akinci p. 3)
    ! -------------------------------------------------------------------------------------------------------
    ! polynomes de degré 2 et 3
    function C_new_1(q, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: q, R_SPH

        ! return
        real(rp) :: C

        ! variables locales
        real(rp), save :: cste = -9185000000.0_rp
        real(rp), save :: alpha = 0.33333_rp

        if ((0.0_rp <= q) .and. (q <= 1.0_rp)) then
            C = cste * (q - alpha * R_SPH) * (q - R_SPH)
        else
            C = 0.0_rp
        end if
    end function

    function C_new_2(q, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: q, R_SPH

        ! return
        real(rp) :: C

        ! variables locales
        real(rp), save :: cste = -9185000000.0_rp
        real(rp), save :: alpha = -0.3333333333_rp
        real(rp), save :: beta = 0.3333333333_rp

        if ((0.0_rp <= q) .and. (q <= 1.0_rp)) then
            C = cste * (q - alpha * R_SPH) * (q - beta * R_SPH) * (q - R_SPH)
        else
            C = 0.0_rp
        end if
    end function

    function C_new_3(q, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: q, R_SPH

        ! return
        real(rp) :: C

        !C = 0.67_rp*q**3 - 1.51_rp*q**2 + 1.03_rp*q - 0.18_rp
        C = 0.6666665624118165_rp*q**3 - 1.513333124823633_rp*q**2 + 1.0266665624118165_rp*q - 0.18_rp
    end function

    ! tests ressemblance fonction akinci
    function C_new_4(q, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: q, R_SPH

        ! return
        real(rp) :: C

        ! variables locales

        C = 10.94_rp*q**4 - 21.17_rp*q**3 + 11.52_rp*q**2 - 1.12_rp*q - 0.18_rp
    end function

    !p(0) = -0.18
    !p(0.27) = 0
    !p(0.5) = 0.18
    !p(1) = 0
    !p'(0) = 0
    !p'(0.5) = 0
    !p'(1) = 0
    function C_new_5(q, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: q, R_SPH

        ! return
        real(rp) :: C

        ! variables locales

        C = -28.56578847417548_rp*q**6 + 81.37736542252645_rp*q**5 - 77.71881254107033_rp*q**4 + &
            25.208682711263233_rp*q**3 - 0.12144711854387105_rp*q**2 + 0.0_rp*q - 0.18_rp
    end function

    function C_new_6(q, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: q, R_SPH

        ! return
        real(rp) :: C

        if ((0.0_rp <= q) .and. (q <= 1.0_rp/sqrt(2.0_rp))) then
            C = sqrt(q) -0.18_rp
        else if ((1.0_rp/sqrt(2.0_rp) <= q) .and. (q <= 1.0_rp)) then
            C = -6.99313528071068_rp*q**5 + 174.83919617235026_rp*q**4 - 481.42706065296426_rp*q**3 &
                + 504.05402149427505_rp*q**2 - 228.2179693155049_rp*q + 37.74494758255456_rp
        end if
    end function

    ! tests spheric
    function C_new_7(q, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: q, R_SPH

        ! return
        real(rp) :: C

        ! variables locales
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Force de cohésion des particules i <- j
    ! -------------------------------------------------------------------------------------------------------
    subroutine F_ij_new(sigma, i, j, part, F)
        ! paramètres
        real(rp), intent(in) :: sigma
        integer, intent(in) :: i, j
        type(Particules), intent(in) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
            ! cohesion force
            !TODO : vérifier si il y a un moins devant le sigma
            F = -sigma * part%w(i) * part%R * C_new_5(fnorme2(part%x(i, :) - part%x(j, :)) / part%R, part%R) &
                * (part%x(i, :) - part%x(j, :)) / fnorme2(part%x(i, :) - part%x(j, :))

            ! curvature force
            !F = F - sigma * part%R * (part%gradR(i, :) - part%gradR(j, :))
        else
            F = 0.0_rp
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Force de tension de surface solver SPH
    ! -------------------------------------------------------------------------------------------------------
    ! sigma : coefficient de tension de surface (gamma)
    ! i : numéro d'une particule
    ! part : liste des particules
    subroutine FTS_new(sigma, i, part, F)
        ! paramètres
        real(rp), intent(in) :: sigma
        integer, intent(in) :: i
        type(Particules), intent(inout) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        ! variables locales
        real(rp), dimension(SPH_D) :: ni, nj, Fij
        integer :: k, j

        F = 0.0_rp
        do j = 1, i - 1
            call F_ij_new(sigma, i, j, part, Fij)
            F = F + Fij
        end do
        do j = i + 1, part%n
            call F_ij_new(sigma, i, j, part, Fij)
            F = F + Fij
        end do
    end subroutine



    ! =======================================================================================================
    ! FTS RAYON
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! Force de cohésion des particules i <- j
    ! -------------------------------------------------------------------------------------------------------
    subroutine F_ij_rayon(sigma, i, j, part, F)
        ! paramètres
        real(rp), intent(in) :: sigma
        integer, intent(in) :: i, j
        type(Particules), intent(in) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        ! variables locales
        integer :: k
        real(rp), dimension(SPH_D) :: grad
        real(rp) :: prod, div_nori

        if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
            F = -4.0_rp * sigma * part%div_nor(i) * part%w(j) * part%dWij(i, j, :)
        else
            F = 0.0_rp
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! div(ni) = div(gradR)
    ! -------------------------------------------------------------------------------------------------------
    subroutine set_div_nor(part, i)
        ! paramètres
        type(Particules), intent(inout) :: part
        integer, intent(in) :: i

        ! variables locales
        integer :: k
        real(rp) :: prod

        if (fnorme2(part%x(i, :) - part%centre) <= part%rayon - part%R + part%dx/2.0_rp) then
            ! particules intérieurs
            part%div_nor(i) = 0.0_rp
        else
            ! particules périphériques
            part%div_nor(i) = 0.0_rp
            do k = 1, part%n
                call prodScal(part%nor(k, :) - part%nor(i, :), part%dWij(i, k, :), prod)
                part%div_nor(i) = part%div_nor(i) + part%w(k) * prod
            end do
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Force de tension de surface solver SPH
    ! -------------------------------------------------------------------------------------------------------
    ! sigma : coefficient de tension de surface (gamma)
    ! i : numéro d'une particule
    ! part : liste des particules
    subroutine FTS_rayon(sigma, i, part, F)
        ! paramètres
        real(rp), intent(in) :: sigma
        integer, intent(in) :: i
        type(Particules), intent(inout) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        ! variables locales
        real(rp), dimension(SPH_D) :: ni, nj, Fij
        integer :: k, j

        call set_div_nor(part, i)
        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if (part%div_nor(i) /= 0.0_rp) then
            write (*, '("# i :",1I4,"     # div(ni) :",1F9.3,"     # 1/Rg :",1F6.1)'), i, -part%div_nor(i), 1.0_rp / 0.01_rp
        end if
        ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        F = 0.0_rp
        do j = 1, i - 1
            call F_ij_rayon(sigma, i, j, part, Fij)
            F = F + Fij
        end do
        do j = i + 1, part%n
            call F_ij_new(sigma, i, j, part, Fij)
            F = F + Fij
        end do
    end subroutine



    ! =======================================================================================================
    ! CALCULS COMMUNS TENSION DE SURFACE
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! calculer la tension de surface pour toute particule i
    ! -------------------------------------------------------------------------------------------------------
    ! sigma : coefficient de tension superficielle (cste)
    ! part : liste des particules
    subroutine set_fts(FTS_func, sigma, part)
        ! paramètres
        interface
            subroutine FTS_func(sigma_, i_, part_, F)
                use var
                real(rp), intent(in) :: sigma_
                integer, intent(in) :: i_
                type(Particules), intent(inout) :: part_
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

END MODULE tsup
