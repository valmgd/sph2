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

    character(len=5), parameter :: scenario_a = "pave "
    character(len=5), parameter :: scenario_b = "bulle"
    character(len=5), parameter :: choix = scenario_b

    ! nombre de points par axes, indice de boucle
    integer :: n, i
    ! bornes du domaine dans chaque direction
    real(rp), dimension(d, 2) :: bornes
    ! pour une bulle
    real(rp), dimension(d) :: centre
    real(rp) :: rayon

    ! particules
    type(Particules) :: p



    ! lecture du fichier d'entrée
    call readValues("../entrees/constantes", n, bornes, centre, rayon)

    ! création du maillage initial
    select case (choix)
    case (scenario_a)
        ! maillage d'un carré / cube
        call pave(d, n, bornes, "../sorties/x", p)
        call init_pression_pave(bornes, p)
    case (scenario_b)
        ! maillage d'une bulle
        call bulle(d, n, centre, rayon, "../sorties/x", p)
        call init_pression_bulle(centre, rayon, p)
    case default
        write (*, *) "choix de scénario invalide"
    end select

    call normale_surface_GR(p, "../sorties/normale.dat")



    ! *******************************************************************************************************
    call rm_Particules(p)

END PROGRAM main
