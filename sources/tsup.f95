! ===========================================================================================================
! Module pour tension superficielle
!
!   function C_Akinci(r, R_SPH) result(C)
!   subroutine FTS_akinci(sigma, i, part, F)
!   subroutine F_pp_kordilla(y, R_SPH, i, j, x, F)
!   subroutine F_TS_kordilla()
!   subroutine F_TS_cohesion(y, i, part, F)
!   subroutine set_fts(FTS_func, sigma, part)
! ===========================================================================================================

MODULE tsup

    use var

    implicit none

contains

    ! =======================================================================================================
    ! TENSION SUPERFICIELLE AKINCI
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



    ! =======================================================================================================
    ! TENSION SUPERFICIELLE KORDILLA
    ! =======================================================================================================

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
                use math
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