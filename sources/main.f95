! ===========================================================================================================
! Code SPH 2D
!
!
! ===========================================================================================================

PROGRAM main

    use math
    use var
    use sph
    use donnees

    implicit none

    ! choix de scénario
    character(len=*), parameter :: scenario_a = "pave"
    character(len=*), parameter :: scenario_b = "bulle"
    character(len=*), parameter :: choix = scenario_b

    ! nombre de points par axes, indice de boucle
    integer :: n, i, d_Omega
    ! bornes du domaine dans chaque direction
    real(rp), dimension(:, :), allocatable :: bornes
    ! pour une bulle
    real(rp), dimension(:), allocatable :: centre
    real(rp) :: rayon, intervalle, sigma

    ! particules
    type(Particules) :: p



    ! lecture du fichier d'entrée
    call readValues("../entrees/constantes", d_Omega, sigma, intervalle, n, bornes, centre, rayon)

    ! création du maillage initial
    select case (choix)
    case (scenario_a)
        ! maillage d'un carré / cube
        call pave(SPH_D, n, bornes, "../sorties/x", p)
        call init_var_pave(bornes, p)
    case (scenario_b)
        ! maillage d'une bulle
        call bulle(SPH_D, n, centre, rayon, "../sorties/x", p)
        call init_var_bulle(centre, rayon, p)
        print *, maxval(p%P)
    case default
        write (*, *) "choix de scénario invalide"
    end select

    call normale_surface_GR(p, "../sorties/normale.dat")
    call write_var(p, "../sorties/rho.dat", "../sorties/u.dat", "../sorties/P.dat")



    ! tension de surface
    call set_gradR(p)
    call set_fts(FTS_akinci, DONNEES_SIGMA, p)

    ! schéma SPH (équation 2)
    call iter_SPH(p)



    ! *******************************************************************************************************
    call system("gnuplot ../sources/Plot.gnu")
    call rm_Particules(p)
    deallocate(bornes, centre)

END PROGRAM main
