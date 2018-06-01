! ===========================================================================================================
! Noyaux SPH
! Définition de différent noyau SPH et de leur gradient. Contient deux pointeurs sur fonctions, globaux au
! module sph, et un "seter" pour choisir le noyau voulu au début du code.
!
!   subroutine set_W_SPH(choix)
!   function theta(q)
!   function W_SPH_wendland(x, R) result(W)
!   subroutine dx_W_SPH_wendland(z, R, grad)
!   function W_SPH_kordilla(x, R) result(W)
!   subroutine dx_W_SPH_kordilla(z, R, grad)
! ===========================================================================================================

MODULE kernel

    use var

    implicit none

    abstract  interface
        function W_SPH_version(x, R) result(W)
            use var
            real(rp), dimension(SPH_D), intent(in) :: x
            real(rp), intent(in) :: R
            real(rp) :: W
        end function
    end interface

    procedure(W_SPH_version), pointer :: W_SPH => NULL()


    abstract interface
        subroutine dx_W_SPH_version(x, R, grad)
            use var
            real(rp), dimension(SPH_D), intent(in) :: x
            real(rp), intent(in) :: R
            real(rp), dimension(SPH_D), intent(out) :: grad
        end subroutine
    end  interface

    procedure(dx_W_SPH_version), pointer :: dx_W_SPH => NULL()

contains

    ! -------------------------------------------------------------------------------------------------------
    ! Choix du noyau SPH
    ! -------------------------------------------------------------------------------------------------------
    ! choix : chaîne de caractère parmi :
    !         - "wendland"
    !         - "tartakovsky"
    !         - "liu"
    subroutine set_W_SPH(choix)
        ! paramètres
        character(len=*), intent(in) :: choix

        if (trim(choix) == "wendland") then
            W_SPH => W_SPH_wendland
            dx_W_SPH => dx_W_SPH_wendland
        else if (trim(choix) == "tartakovsky") then
            W_SPH => W_SPH_tartakovsky
            dx_W_SPH => dx_W_SPH_tartakovsky
        else if (trim(choix) == "liu") then
            w_SPH => W_SPH_liu
            dx_W_SPH => dx_W_SPH_liu
        end if
    end subroutine



    ! =======================================================================================================
    ! NOYAU SPH WENDLAND {{{
    ! =======================================================================================================

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
            theta = SPH_NORM_NOYAU_WENDLAND * (1.0_rp - q)**4 * (1.0_rp + 4.0_rp*q)
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
    function W_SPH_wendland(x, R) result(W)
        ! paramètres
        real(rp), dimension(SPH_D), intent(in) :: x
        real(rp), intent(in) :: R

        ! return
        real(rp) :: W

        ! variables locales
        real(rp) :: norm

        call norme2(x, norm)
        if (norm <= R) then
            W = theta(norm / R) / (SPH_VOLUME_SUPP_WENDLAND * R**SPH_D)
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
    subroutine dx_W_SPH_wendland(z, R, grad)
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
            grad = (SPH_NORM_NOYAU_WENDLAND / (SPH_VOLUME_SUPP_WENDLAND * R**SPH_D)) * (-20.0_rp / R) * (z / norm) * q * (1-q)**3
        else
            grad = 0.0_rp
        end if
    end subroutine
    ! }}}



    ! =======================================================================================================
    ! NOYAU SPH TARTAKOVSKY (apparaît dans le papier de Kordilla) {{{
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! noyau Tartakovsky
    ! -------------------------------------------------------------------------------------------------------
    function W_SPH_tartakovsky(x, R) result(W)
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

        else if ((R / 3.0_rp <= length) .and. (length < 2.0_rp * R / 3.0_rp)) then
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
    ! Gradient en x de W(x, y, R) = W(x - y, R)
    ! -------------------------------------------------------------------------------------------------------
    ! z : élt de R**2 (en pratique coordonnées d'une particule)
    ! R : rayon SPH
    ! grad : résultat (R**2)
    subroutine dx_W_SPH_tartakovsky(z, R, grad)
        ! paramètres
        real(rp), dimension(SPH_D), intent(in) :: z
        real(rp), intent(in) :: R
        real(rp), dimension(SPH_D), intent(out) :: grad

        ! variables locales
    end subroutine
    ! }}}



    ! =======================================================================================================
    ! NOYAU SPH LIU (méthode de TS de Kordilla basé dessus) {{{
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! noyau Liu
    ! -------------------------------------------------------------------------------------------------------
    function W_SPH_liu(x, R) result(W)
        ! paramètres
        real(rp), dimension(SPH_D), intent(in) :: x
        real(rp), intent(in) :: R

        ! return
        real(rp) :: W

        ! variables locales
        real(rp) :: length

        length = fnorme2(x)

        if ((0.0_rp <= length) .and. (length < R / 2.0_rp)) then
            W = 1.0_rp - 6.0_rp * (length/R)**2 + 6.0_rp * (length/R)**3
        else if ((R / 2.0_rp <= length) .and. (length < R)) then
            W = 2.0_rp * (1.0_rp - length/R)**3
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
    subroutine dx_W_SPH_liu(z, R, grad)
        ! paramètres
        real(rp), dimension(SPH_D), intent(in) :: z
        real(rp), intent(in) :: R
        real(rp), dimension(SPH_D), intent(out) :: grad

        ! variables locales
    end subroutine
    ! }}}

END MODULE kernel
