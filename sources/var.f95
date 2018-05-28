! ===========================================================================================================
! Module contenant les variables nécessaire au bloc SPH
!
!   subroutine set_SPH_D(valeur)
!   subroutine set_SPH_I(valeur)
!   subroutine rm_Particules(part)
! ===========================================================================================================

MODULE var

    use math

    implicit none

    ! dimension du domaine
    integer, save :: SPH_D
    ! intervalle : R_SPH = SPH_I * dx
    real(rp), save :: SPH_I

    ! constantes de normalisation pour noyaux SPH
    real(rp), save :: SPH_NORM_NOYAU_WENDLAND
    real(rp), save :: SPH_VOLUME_SUPP_WENDLAND

    real(rp), save :: SPH_NORM_NOYAU_TARTAKOVSKY

    ! constante de normalisation pour noyau cohésion Akinci (tension de surface)
    real(rp), save :: SPH_NORM_COHESION_AKINCI



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

        ! noyau
        real(rp), dimension(:, :, :), allocatable :: dWij

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
        SPH_NORM_NOYAU_WENDLAND    =                  7.0_rp**(3-SPH_D) *               14.0_rp**(SPH_D-2)
        SPH_VOLUME_SUPP_WENDLAND   =                      pi**(3-SPH_D) *    (4.0_rp*pi/3.0_rp)**(SPH_D-2)

        SPH_NORM_NOYAU_TARTAKOVSKY =   (63.0_rp/478.0_rp/pi)**(3-SPH_D) * (81.0_rp/359.0_rp/pi)**(SPH_D-2)

        SPH_NORM_COHESION_AKINCI   = (6672.0_rp/185.0_rp/pi)**(3-SPH_D) *          (32.0_rp/pi)**(SPH_D-2)
    end subroutine

    subroutine set_SPH_I(valeur)
        ! paramètres
        real(rp), intent(in) :: valeur
        SPH_I = valeur
    end subroutine



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

        if (allocated(part%dWij)) then
            deallocate(part%dWij)
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

END MODULE var
