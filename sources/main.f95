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

    character(len=4), parameter :: scenario_a = "cube"
    character(len=4), parameter :: scenario_b = "rond"
    character(len=4), parameter :: choix = scenario_a

    ! nombre de points par axes, indice de boucle
    integer :: n, i
    ! bornes du domaine dans chaque direction
    real(rp), dimension(d, 2) :: bornes
    ! pour une bulle
    real(rp), dimension(d) :: centre
    real(rp) :: rayon

    ! particules
    type(Particules) :: p

    ! lecture paramètres
    open(unit = 10, file = "../entrees/constantes")
    read (10, *) n
    read (10, *) centre
    read (10, *) rayon
    do i = 1, d
        read (10, *) bornes(i, :)
    end do
    close(10)

    ! maillage d'un carré / cube
    call cube(d, n, bornes, "../sorties/x", p)

    ! *******************************************************************************************************
    call rm_Particules(p)

END PROGRAM main
