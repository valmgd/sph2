! ===========================================================================================================
! Code SPH 2D / 3D
!
!
! ===========================================================================================================

PROGRAM main

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
    real(rp), dimension(:, :), allocatable :: normales
    real(rp), dimension(:), allocatable :: div_normales
    real(rp) :: angle, prod, rtemp
    integer :: k



    ! =======================================================================================================
    ! Computation
    ! =======================================================================================================
    call cpu_time(t_start)

    ! lecture du fichier d'entrée
    call readValues("../entrees/constantes", d_Omega, sigma, intervalle, n, bornes, centre, rayon)

    ! choix du noyau SPH
    call set_W_SPH("wendland")
    ! call set_W_SPH("liu")

    ! -------------------------------------------------------------------------------------------------------
    ! création du maillage initial
    ! -------------------------------------------------------------------------------------------------------
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

    ! -------------------------------------------------------------------------------------------------------
    ! tension de surface
    ! -------------------------------------------------------------------------------------------------------
    print *, "Set surface tension. Thinking..."
    ! call set_fts(FTS_akinci, DONNEES_SIGMA, p)
    ! call set_fts(FTS_liu, DONNEES_SIGMA, p)
    ! call set_fts(FTS_new, DONNEES_SIGMA, p)
    call set_fts(FTS_rayon, DONNEES_SIGMA, p)



    ! schéma SPH (équation 2)
    print *, "Compute mvt quantity. Thinking..."
    call iter_SPH(p, centre)

    call cpu_time(t_end)




    ! =======================================================================================================
    ! Post traitement : sauvegardes, affichages et vérifs
    ! =======================================================================================================
    call write_Particules(p, '../sorties')

    print *, "________________________________________________"
    print *, "                     |"
    print *, "Nombre de particules |", p%n
    print *, "dx                   |", p%dx
    print *, "Rayon SPH            |", p%R
    print *, "_____________________|__________________________"

    call quarter(p, centre)





    ! Vérifs et sorties de graphes {{{
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! vérif intégrale sur l'intersection des deux disques {{{
    ! centre1 = 0.0_rp
    ! centre2 = p%x(p%n, :)
    ! somme = 0.0_rp
    ! do i = 1, p%n - 1
    !     x = p%x(p%n, :)
    !     y = p%x(i, :)
    !     point = x - y
    !     if ((fnorme2(y - centre1) <= rayon) .and. (fnorme2(y - centre2) <= p%R)) then
    !         call dx_W_SPH(point, p%R, grad)
    !         somme = somme + grad * p%dx**2
    !         ! print *, i, somme, grad
    !     end if
    ! end do
    ! print *, "Intégrale sur l'intersection :", somme
    ! }}}
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! affichage fonction cohésion {{{
    ! p%R = 1.0_rp
    ! rr = linspace(0.0_rp, 1.0_rp, 1000)

    ! akinci
    ! do i = 1, 1000
    !     Ck(i) = C_akinci(rr(i), p%R)
    ! end do
    ! call saveSol(rr, Ck, "../sorties/Cakinci.dat")

    ! liu (kordilla)
    ! do i = 1, 1000
    !     Ck(i) = C_liu(rr(i), p%R)
    !     AW1(i) =  2.0_rp * W_SPH_liu((/ rr(i), 0.0_rp /), 0.8_rp * p%R)
    !     BW2(i) = -1.0_rp * W_SPH_liu((/ rr(i), 0.0_rp /), 1.0_rp * p%R)
    ! end do
    ! call saveSol(rr, AW1, "../sorties/AW1.dat")
    ! call saveSol(rr, BW2, "../sorties/BW2.dat")
    ! call saveSol(rr, Ck, "../sorties/Cliu.dat")

    ! new
    ! do i = 1, 1000
    !     Ck(i) = C_new_5(rr(i), p%R)
    ! end do
    ! call saveSol(rr, Ck, "../sorties/Cnew.dat")
    ! }}}
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! }}}

    ! normales exactes {{{
    allocate(normales(p%n, 2), div_normales(p%n))
    do i = 1, p%n
        if (p%x(i, 1) == 0.0_rp) then
            if (p%x(i, 2) > 0.0_rp) then
                angle = pi / 2.0_rp
            else
                angle = -pi / 2.0_rp
            endif
        else
            angle = atan(p%x(i, 2) / p%x(i, 1))
        endif

        if (p%x(i, 1) >= 0.0_rp) then
            normales(i, :) = (/ -cos(angle), -sin(angle) /)
        else
            normales(i, :) = (/ cos(angle), sin(angle) /)
        endif
    end do
    open(unit = 10, file = '../sorties/normales.dat')
        do i = 1, p%n
            write (10, *) normales(i, :)
        end do
    close(10)

    do i = 1, p%n
        div_normales(i) = 0.0_rp
        do k = 1, p%n
            call prodScal(normales(k, :) - normales(i, :), p%dWij(i, k, :), prod)
            div_normales(i) = div_normales(i) + p%w(k) * prod
        end do
    end do

    do i = 1, p%n
        if (fnorme2(p%x(i, :)) < rayon - p%R) then
            print *, 1.0_rp / fnorme2(p%x(i, :)), -div_normales(i)
        else
            print *, "                                            ", 1.0_rp / fnorme2(p%x(i, :)), -div_normales(i)
        endif
    end do
    ! }}}




    ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    print *
    write (*, '("### temps d''éxecution :",1F6.2)'), t_end - t_start
    call rm_Particules(p)
    deallocate(bornes, centre)
    ! call system('../pyplot/main.py')

END PROGRAM main
