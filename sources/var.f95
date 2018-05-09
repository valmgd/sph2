! ===========================================================================================================
! Module contenant les variables nécessaire au bloc SPH
!
! ===========================================================================================================

MODULE var

    use math

    implicit none

    ! dimension du domaine
    integer, save :: SPH_D
    ! intervalle : R_SPH = SPH_I * dx
    real(rp), save :: SPH_I

    ! constantes pour noyau SPH
    real(rp), save :: SPH_C_NORM_NOYAU
    real(rp), save :: SPH_VOLUME_SUPP

    ! constantes pour fonction Akinci (tension de surface)
    real(rp), save :: SPH_C_NORM_AKINCI


    type :: Particules
        ! nombre de particules
        integer :: n
        ! coordonnées des particules
        real(rp), dimension(:, :), allocatable :: x
        ! volume des particules
        real(rp), dimension(:), allocatable :: w

        ! rayon SPH
        real(rp) :: R
        real(rp) :: dx

        ! gradient de la fc cste R
        real(rp), dimension(:, :), allocatable :: gradR
        real(rp), dimension(:, :), allocatable :: fts

        ! variables équations d'Euler
        real(rp), dimension(:), allocatable :: rho
        real(rp), dimension(:, :), allocatable :: u
        real(rp), dimension(:), allocatable :: P
    end type Particules

contains

    subroutine set_SPH_D(valeur)
        ! paramètres
        integer, intent(in) :: valeur
        SPH_D = valeur
        SPH_C_NORM_NOYAU  =              7.0_rp**(3-SPH_D) *            14.0_rp**(SPH_D-2)
        SPH_VOLUME_SUPP   =                  pi**(3-SPH_D) * (4.0_rp*pi/3.0_rp)**(SPH_D-2)
        SPH_C_NORM_AKINCI = (417.0_rp/370.0_rp)**(3-SPH_D) *             1.0_rp**(SPH_D-2)
    end subroutine

    subroutine set_SPH_I(valeur)
        ! paramètres
        real(rp), intent(in) :: valeur
        SPH_I = valeur
    end subroutine

END MODULE var
