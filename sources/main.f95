! ===========================================================================================================
! Code SPH 2D
!
!
! ===========================================================================================================

PROGRAM main

    use math
    use donnees
    use sph

    implicit none

    character(len=4), parameter :: choix_a = "cube"
    character(len=4), parameter :: choix_b = "rond"
    character(len=4) :: choix = choix_a

    ! nombre de points par axes
    integer :: n
    ! bornes du domaine dans chaque direction
    real(rp) :: xmin, xmax
    ! liste des particules
    real(rp), dimension(:, :), allocatable :: x
    ! volume des particules
    real(rp), dimension(:), allocatable :: w

    ! lecture paramètres
    open(unit = 10, file = "../entrees/constantes")
    read (10, *) n
    read (10, *) xmin
    read (10, *) xmax
    close(10)

    ! maillage d'un carré / cube
    call cube(d, n, xmin, xmax, "../sorties/x", x, w)

    ! *******************************************************************************************************
    deallocate(x)

END PROGRAM main
