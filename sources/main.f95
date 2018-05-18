! ===========================================================================================================
! Code SPH 2D / 3D
!
!
! ===========================================================================================================

PROGRAM main

    ! libraries
    use sph

    ! modules perso
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

    real(rp), dimension(1000) :: rr, Ck, AW1, BW2



    ! lecture du fichier d'entrée
    call readValues("../entrees/constantes", d_Omega, sigma, intervalle, n, bornes, centre, rayon)

    ! choix du noyau SPH
    call set_W_SPH("wendland")
    !call set_W_SPH("liu")

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
    case default
        write (*, *) "choix de scénario invalide"
    end select

    call normale_surface_GR(p, "../sorties/normale.dat")
    call write_var(p, "../sorties/rho.dat", "../sorties/u.dat", "../sorties/P.dat")



    ! tension de surface
    !call set_fts(FTS_akinci, DONNEES_SIGMA, p)
    !call set_fts(FTS_liu, DONNEES_SIGMA, p)
    call set_fts(FTS_new, DONNEES_SIGMA, p)

    ! schéma SPH (équation 2)
    call iter_SPH(p, centre)
    call quarter(p, centre)
    print *, "np", p%n


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print *, "# dx", p%dx
    p%R = 1.0_rp
    rr = linspace(0.0_rp, 1.0_rp, 1000)
    do i = 1, 1000
        Ck(i) = C_liu(rr(i), p%R)
        AW1(i) =  2.0_rp * W_SPH_liu((/ rr(i), 0.0_rp /), 0.8_rp * p%R)
        BW2(i) = -1.0_rp * W_SPH_liu((/ rr(i), 0.0_rp /), 1.0_rp * p%R)
    end do
    call saveSol(rr, Ck, "../sorties/Ckordilla.dat")
    call saveSol(rr, AW1, "../sorties/AW1.dat")
    call saveSol(rr, BW2, "../sorties/BW2.dat")
    print *, "# R", p%R

    do i = 1, 1000
        Ck(i) = C_akinci(rr(i), p%R)
    end do
    call saveSol(rr, Ck, "../sorties/Cakinci.dat")

    do i = 1, 1000
        Ck(i) = C_new_2(rr(i), p%R)
    end do
    call saveSol(rr, Ck, "../sorties/Cnew.dat")



    ! *******************************************************************************************************
    call system("gnuplot ../sources/Plot.gnu")
    call rm_Particules(p)
    deallocate(bornes, centre)

END PROGRAM main
