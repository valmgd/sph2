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

        ! rayon SPH et pas d'espace maillage initial
        real(rp) :: R
        real(rp) :: dx

        ! coordonnées et volume des particules
        real(rp), dimension(:, :), allocatable :: x
        real(rp), dimension(:), allocatable :: w

        ! variables équations d'Euler
        real(rp), dimension(:), allocatable :: rho
        real(rp), dimension(:, :), allocatable :: u
        real(rp), dimension(:), allocatable :: P

        ! gradient de la fc cste R, normale à la surface, divergence de la normale
        real(rp), dimension(:, :), allocatable :: gradR
        real(rp), dimension(:, :), allocatable :: nor
        real(rp), dimension(:), allocatable :: div_nor

        ! gradient de pression, tension superficielle, quantité de mouvement
        real(rp), dimension(:, :), allocatable :: grad_P
        real(rp), dimension(:, :), allocatable :: fts
        real(rp), dimension(:, :), allocatable :: dmu

        ! noyau
        real(rp), dimension(:, :, :), allocatable :: dWij

        ! infos surface libre
        real(rp), dimension(:), allocatable :: centre
        real(rp) :: rayon
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

        ! coordonnées et volume des particules
        if (allocated(part%x)) then
            deallocate(part%x)
        end if
        if (allocated(part%w)) then
            deallocate(part%w)
        end if

        ! variables équations d'Euler
        if (allocated(part%rho)) then
            deallocate(part%rho)
        end if
        if (allocated(part%u)) then
            deallocate(part%u)
        end if
        if (allocated(part%P)) then
            deallocate(part%P)
        end if

        ! gradient de la fc cste R, normale à la surface, divergence de la normale
        if (allocated(part%gradR)) then
            deallocate(part%gradR)
        end if
        if (allocated(part%nor)) then
            deallocate(part%nor)
        end if
        if (allocated(part%div_nor)) then
            deallocate(part%div_nor)
        end if

        ! gradient de pression, tension superficielle, quantité de mouvement
        if (allocated(part%grad_P)) then
            deallocate(part%grad_P)
        end if
        if (allocated(part%fts)) then
            deallocate(part%fts)
        end if
        if (allocated(part%dmu)) then
            deallocate(part%dmu)
        end if

        ! noyau
        if (allocated(part%dWij)) then
            deallocate(part%dWij)
        end if

        ! infos surface libre
        if (allocated(part%centre)) then
            deallocate(part%centre)
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! -------------------------------------------------------------------------------------------------------
    subroutine write_Particules(part, dir)
        ! paramètres
        type(Particules), intent(in) :: part
        character(len=*), intent(in) :: dir

        ! variables locales
        integer :: i

        open(unit = 10, file = dir // '/x.dat')
        open(unit = 20, file = dir // '/w.dat')
        open(unit = 30, file = dir // '/rho.dat')
        open(unit = 40, file = dir // '/u.dat')
        open(unit = 50, file = dir // '/P.dat')
        open(unit = 60, file = dir // '/nor.dat')
        open(unit = 70, file = dir // '/div_nor.dat')
        open(unit = 80, file = dir // '/grad_P.dat')
        open(unit = 90, file = dir // '/fts.dat')
        open(unit = 100, file = dir // '/dmu.dat')
        do i = 1, part%n
            write (10, *) part%x(i, :)
            write (20, *) part%w(i)

            write (30, *) part%rho(i)
            write (40, *) part%u(i, :)
            write (50, *) part%P(i)

            write (60, *) part%nor(i, :)
            write (70, *) part%div_nor(i)

            write (80, *) part%grad_P(i, :)
            write (90, *) part%fts(i, :)
            write (100, *) part%dmu(i, :)
        end do
        close(10)
        close(20)
        close(30)
        close(40)
        close(50)
        close(60)
        close(70)
        close(80)
        close(90)
        close(100)
    end subroutine

END MODULE var
