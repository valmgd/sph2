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
    ! TESTS DIVERS
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! Noyau cohésion pour tension de surface (cf C(r) Akinci p. 3)
    ! -------------------------------------------------------------------------------------------------------
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
            F = -sigma * part%w(i) * part%R * C_new_2(fnorme2(part%x(i, :) - part%x(j, :)) / part%R, part%R) &
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
