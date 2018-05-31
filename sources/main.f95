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
    real(rp) :: rayon, intervalle, sigma, t_start, t_end

    ! particules
    type(Particules) :: p

    !!!!!!!!!!!!!!!!!
    real(rp), dimension(1000) :: rr, Ck, AW1, BW2
    real(rp), dimension(2) :: grad, point, centre1, centre2, x, y, somme



    ! -------------------------------------------------------------------------------------------------------
    ! Computation
    ! -------------------------------------------------------------------------------------------------------
    call cpu_time(t_start)

    ! lecture du fichier d'entrée
    call readValues("../entrees/constantes", d_Omega, sigma, intervalle, n, bornes, centre, rayon)

    ! choix du noyau SPH
    call set_W_SPH("wendland")
    !call set_W_SPH("liu")

    ! création du maillage initial
    print *, "Create mesh. Thinking..."
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
    print *, "Set surface tension. Thinking..."
    call set_fts(FTS_akinci, DONNEES_SIGMA, p)
    ! call set_fts(FTS_liu, DONNEES_SIGMA, p)
    ! call set_fts(FTS_new, DONNEES_SIGMA, p)
    ! call set_fts(FTS_rayon, DONNEES_SIGMA, p)

    ! schéma SPH (équation 2)
    print *, "Compute mvt quantity. Thinking..."
    call iter_SPH(p, centre)

    call cpu_time(t_end)



    ! -------------------------------------------------------------------------------------------------------
    ! Affichages et vérifs
    ! -------------------------------------------------------------------------------------------------------
    print *, "________________________________________________"
    print *, "                     |"
    print *, "Nombre de particules |", p%n
    print *, "dx                   |", p%dx
    print *, "Rayon SPH            |", p%R
    print *, "_____________________|__________________________"
    call quarter(p, centre)










    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! vérif intégrale sur l'intersection des deux disques
    centre1 = 0.0_rp
    centre2 = p%x(p%n, :)
    somme = 0.0_rp
    do i = 1, p%n - 1
        x = p%x(p%n, :)
        y = p%x(i, :)
        point = x - y
        if ((fnorme2(y - centre1) <= rayon) .and. (fnorme2(y - centre2) <= p%R)) then
            call dx_W_SPH(point, p%R, grad)
            somme = somme + grad * p%dx**2
            ! print *, i, somme, grad
        end if
    end do
    print *, "Intégrale sur l'intersection :", somme
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! affichage fonction cohésion
    p%R = 1.0_rp
    rr = linspace(0.0_rp, 1.0_rp, 1000)

    ! akinci
    do i = 1, 1000
        Ck(i) = C_akinci(rr(i), p%R)
    end do
    call saveSol(rr, Ck, "../sorties/Cakinci.dat")

    ! liu (kordilla)
    do i = 1, 1000
        Ck(i) = C_liu(rr(i), p%R)
        AW1(i) =  2.0_rp * W_SPH_liu((/ rr(i), 0.0_rp /), 0.8_rp * p%R)
        BW2(i) = -1.0_rp * W_SPH_liu((/ rr(i), 0.0_rp /), 1.0_rp * p%R)
    end do
    call saveSol(rr, AW1, "../sorties/AW1.dat")
    call saveSol(rr, BW2, "../sorties/BW2.dat")
    call saveSol(rr, Ck, "../sorties/Cliu.dat")

    ! new
    do i = 1, 1000
        Ck(i) = C_new_5(rr(i), p%R)
    end do
    call saveSol(rr, Ck, "../sorties/Cnew.dat")
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





    ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    print *
    write (*, '("### temps d''éxecution :",1F6.2)'), t_end - t_start
    call rm_Particules(p)
    deallocate(bornes, centre)
    write (*, '(/,/,/,"_______________",/,"CALLING GNUPLOT")')
    call system("gnuplot ../sources/Plot.gnu")
    write (*, '("Done.")')

END PROGRAM main
