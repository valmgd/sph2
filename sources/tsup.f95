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
    function C_akinci(z, R_SPH) result(C)
        ! paramètres
        real(rp), intent(in) :: z, R_SPH

        ! return
        real(rp) :: C

        if ((0.0_rp < z) .and. (z <= R_SPH / 2.0_rp)) then
            C = SPH_NORM_COHESION_AKINCI * (32.0_rp / (pi * R_SPH**SPH_D * R_SPH**6)) * &
                (2.0_rp * (R_SPH - z)**3 * z**3 - R_SPH**6 / 64.0_rp)

        else if ((R_SPH / 2.0_rp < z) .and. (z <= R_SPH)) then
            C = SPH_NORM_COHESION_AKINCI * (32.0_rp / (pi * R_SPH**SPH_D * R_SPH**6)) * &
                (R_SPH - z)**3 * z**3

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
            F = sigma * part%w(i) * part%R * C_akinci(fnorme2(part%x(i, :) - part%x(j, :)), part%R) &
                * (part%x(i, :) - part%x(j, :)) / fnorme2(part%x(i, :) - part%x(j, :))

            ! curvature force
            F = F + sigma * part%R * (part%gradR(i, :) - part%gradR(j, :))
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
        real(rp), dimension(SPH_D) :: ni, nj
        integer :: k, j

        F = 0.0_rp
        do j = 1, i - 1
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                F = F + sigma * part%w(i) * part%R * C_akinci(fnorme2(part%x(i, :) - part%x(j, :)), part%R) &
                    * (part%x(i, :) - part%x(j, :)) / fnorme2(part%x(i, :) - part%x(j, :)) &
                    + sigma * part%R * (part%gradR(i, :) - part%gradR(j, :))
            end if
        end do
        do j = i + 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                F = F + sigma * part%w(i) * part%R * C_akinci(fnorme2(part%x(i, :) - part%x(j, :)), part%R) &
                    * (part%x(i, :) - part%x(j, :)) / fnorme2(part%x(i, :) - part%x(j, :)) &
                    + sigma * part%R * (part%gradR(i, :) - part%gradR(j, :))
            end if
        end do
    end subroutine



    ! =======================================================================================================
    ! TENSION SUPERFICIELLE KORDILLA
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! Force d'une particule sur une autre (Kordilla) cohésion ?
    ! -------------------------------------------------------------------------------------------------------
    ! sigma : coeff de tension de surface (gamma)
    !function C_akinci(z, R_SPH) result(C)
    subroutine F_ij_kordilla(sigma, i, j, part, F)
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
            F = sigma * (A * W_SPH_tartakovsky(r, h1) * r / length + B * W_SPH_tartakovsky(r, h2) * r / length)
        else
            F = 0.0_rp
        end if
    end subroutine

    !function C_akinci(z, R_SPH) result(C)
    function C_kordilla(z, R_SPH, sigma) result(C)
        ! paramètres
        real(rp), intent(in) :: z, sigma, R_SPH

        ! return
        real(rp) :: C

        ! variables locales
        real(rp), dimension(SPH_D) :: F
        real(rp) :: s_ff, A, B, h1, h2
        real(rp), dimension(SPH_D) :: r
        real(rp) :: length

        A = 2.0_rp
        B = -1.0_rp
        h1 = 0.8_rp * R_SPH
        h2 = 1.0_rp * R_SPH

        r = (/ z, 0.0_rp /)
        length = fnorme2(r)

        if (length <= R_SPH) then
            C = sigma * (A * W_SPH_tartakovsky(r, h1) * z / z + B * W_SPH_tartakovsky(r, h2) * z / z)
            C = sigma * (A * W_SPH_liu(r, h1) * z / z + B * W_SPH_liu(r, h2) * z / z)
        else
            C = 0.0_rp
        end if

    end function



    ! -------------------------------------------------------------------------------------------------------
    ! force de tension de surface Kordilla
    ! -------------------------------------------------------------------------------------------------------
    subroutine FTS_kordilla(sigma, i, part, F)
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

        do j = 1, part%n
            call F_ij_kordilla(sigma, i, j, part, F)
            F = F + part%w(j) * Fij
        end do
    end subroutine



    ! =======================================================================================================
    ! TESTS DIVERS
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! Force de tension de surface idem akinci mais avec juste F_cohesion
    ! -------------------------------------------------------------------------------------------------------
    ! sigma : coefficient de tension de surface (gamma)
    ! i : numéro d'une particule
    ! part : liste des particules
    subroutine F_TS_cohesion(sigma, i, part, F)
        ! paramètres
        real(rp), intent(in) :: sigma
        integer, intent(in) :: i
        type(Particules), intent(in) :: part
        real(rp), dimension(SPH_D), intent(out) :: F

        ! variables locales
        integer :: k, j

        F = 0.0_rp
        do j = 1, i - 1
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                F = F + sigma * part%w(i) * part%R * C_akinci(fnorme2(part%x(i, :) - part%x(j, :)), part%R) &
                    * (part%x(i, :) - part%x(j, :)) / fnorme2(part%x(i, :) - part%x(j, :))
            end if
        end do
        do j = i + 1, part%n
            if (fnorme2(part%x(i, :) - part%x(j, :)) <= part%R) then
                F = F + sigma * part%w(i) * part%R * C_akinci(fnorme2(part%x(i, :) - part%x(j, :)), part%R) &
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
