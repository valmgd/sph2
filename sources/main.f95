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

    integer, parameter :: scenario = 1

    ! PARTIE 1 : MAILLAGE D'UN PAVÉ 2D
    ! variables pour "maillage" initial
    integer :: n
    real(rp) :: xmin, xmax, ymin, ymax
    real(rp), dimension(:), allocatable :: x1, x2
    real(rp) :: dx, dy

    ! nombre de particules
    integer :: np
    ! coordonnées des particules
    real(rp), dimension(:, :), allocatable :: x
    ! volume des particules
    real(rp), dimension(:), allocatable :: w
    ! rayon SPH
    real(rp) :: R

    ! variables temporaires
    real(rp)               :: real1
    real(rp), dimension(2) :: real2, real2_2, real2_3
    integer, dimension(2) :: sh

    integer :: i, k

    ! vecteurs n de kappa = (1 / Rc) div(n)  (rayon de courbure) normal à la surface libre
    real(rp), dimension(:, :), allocatable :: nvec, grad_P, plot_vec


    ! PARTIE 2 : MAILLAGE D'UNE BULLE
    integer :: cpt
    real(rp) :: rayon
    real(rp), dimension(2) :: centre
    real(rp), dimension(100, 2) :: cc
    real(rp), dimension(100) ::xx

    ! VARIABLES ÉQUATION
    real(rp), dimension(:), allocatable :: P
    real(rp), dimension(:, :), allocatable :: sol_P, nor, fts, d_rwu_dt
    real(rp), dimension(2) :: temp



    open(unit = 10, file = "../entrees/constantes")
    read(10, *) n
    close(10)

    select case (scenario)
    case (0)
        print *, "CARRÉ"
        ! ===================================================================================================
        ! À partir d'un carré
        ! ===================================================================================================

        ! ---------------------------------------------------------------------------------------------------
        ! lecture et initialisation des variables
        ! ---------------------------------------------------------------------------------------------------
        open(unit = 10, file = "../entrees/constantes")
        read(10, *)
        read(10, *) xmin
        read(10, *) xmax
        read(10, *) ymin
        read(10, *) ymax
        close(10)

        ! subdivision des deux dimensions
        x1 = linspace(xmin, xmax, n)
        x2 = linspace(ymin, ymax, n)

        ! création tableau de particules et pas d'espace
        call meshgrid(x1, x2, x)
        call writeMat(x, "../sorties/x.dat")
        dx = x1(2) - x1(1)
        dy = x2(2) - x2(1)

        ! nombre de particules
        sh = shape(x)
        np = sh(1)

        ! ---------------------------------------------------------------------------------------------------
        ! allocations
        ! ---------------------------------------------------------------------------------------------------
        ! vecteur des volumes
        allocate(w(np))
        ! vérif noyau sph
        allocate(nvec(np, 4))
        ! pression
        allocate(P(np))
        ! particules + pression pour écriture fichier
        allocate(sol_P(np, 3))
        ! gradient de pression
        allocate(grad_P(np, 2))
        ! tableau de 4 col pour représentation vecteurs gnuplot
        allocate(plot_vec(np, 4))
        allocate(nor(np, 2))
        ! force de tension de surface
        allocate(fts(np, 2))
        allocate(d_rwu_dt(np, 2))

        w = dx * dy

        ! rayon noyau SPH
        R = 4.0_rp * dx

        xx = linspace(0.00001_rp, R, 100)
        do i = 1, 100
            cc(i, :) = (/ xx(i), C_Akinci(xx(i), R) /)
        end do
        call writeMat(cc, "../sorties/CAkinci.dat")

        ! contour du domaine que nous avons maillé pour graphique
        open(unit = 10, file = "../sorties/cc.dat")
        write (10, *) xmin - 0.5_rp * dx, ymin - 0.5_rp * dy
        write (10, *) xmax + 0.5_rp * dx, ymin - 0.5_rp * dy
        write (10, *) xmax + 0.5_rp * dx, ymax + 0.5_rp * dy
        write (10, *) xmin - 0.5_rp * dx, ymax + 0.5_rp * dy
        write (10, *) xmin - 0.5_rp * dx, ymin - 0.5_rp * dy
        close(10)


        ! ---------------------------------------------------------------------------------------------------
        ! initialisation inconnues équation
        ! ---------------------------------------------------------------------------------------------------
        ! initialisation de la pression
        P = 0.0_rp
        sol_P(:, 1:2) = x
        sol_P(:, 3) = P
        call writeMat(sol_P, "../sorties/P.dat")


        ! approximation du gradient de pression avec l'opérateur GR_p
        do i = 1, np
            call GR_p(i, x, w, R, P, grad_P(i, :))
            ! formattage pour gnuplot
            plot_vec(i, :) = (/ x(i, :), (grad_P(i, :) / fnorme2(grad_P(i, :))) * 0.5_rp * dx /)
        end do
        call writeMat(plot_vec, "../sorties/grad_P.dat")



    case (1)
        print *, "BULLE"
        ! ===================================================================================================
        ! À partir d'une bulle
        ! ===================================================================================================

        ! ---------------------------------------------------------------------------------------------------
        ! maillage de la bulle
        ! ---------------------------------------------------------------------------------------------------
        ! données du disque à mailler
        centre = (/ 0.0_rp, 0.0_rp /)
        rayon = 0.001_rp


        ! maillage du disque et actualisation des paramètre dépendant du maillage
        call meshCircle(centre, rayon, n, x)
        dx = x(2, 1) - x(1, 1)
        R = 4.0_rp * dx
        call writeMat(x, "../sorties/x_circle.dat")
        sh = shape(x)
        np = sh(1)


        ! ---------------------------------------------------------------------------------------------------
        ! allocations
        ! ---------------------------------------------------------------------------------------------------
        ! volume des particules
        allocate(w(np))
        ! vérif noyau sph
        allocate(nvec(np, 4))
        ! pression
        allocate(P(np))
        ! particules + pression pour écriture fichier
        allocate(sol_P(np, 3))
        ! gradient de pression
        allocate(grad_P(np, 2))
        ! tableau de 4 col pour représentation vecteurs gnuplot
        allocate(plot_vec(np, 4))
        allocate(nor(np, 2))
        ! force de tension de surface
        allocate(fts(np, 2))
        allocate(d_rwu_dt(np, 2))

        !w = pi * rayon**2 / real(np, rp)
        w = dx**2


        ! vérification sur l'approx régularisée de la fonction 1 et de son gradient (0, 0)
        i = nint(np / 2.0_rp)
        call AR(i, x, w, R, (/ (1.0_rp, k=1, np) /), real1)
        call GR(i, x, w, R, (/ (1.0_rp, k=1, np) /), real2)
        print *, "vérif AR :", real1
        print *, "vérif GR :", real2


        ! contour de la bulle que nous avons maillé pour graphique
        xx = linspace(0.0_rp, 2.0_rp * pi, 100)
        cc(:, 1) = centre(1) + rayon * cos(xx)
        cc(:, 2) = centre(2) + rayon * sin(xx)
        call writeMat(cc, "../sorties/cc.dat")


        ! ---------------------------------------------------------------------------------------------------
        ! initialisation inconnues équation
        ! ---------------------------------------------------------------------------------------------------
        ! initialisation de la pression
        call init_pression(x, centre, rayon, P)
        sol_P(:, 1:2) = x
        sol_P(:, 3) = P
        call writeMat(sol_P, "../sorties/P.dat")


        ! approximation du gradient de pression avec l'opérateur GR_p
        do i = 1, np
            call GR_p(i, x, w, R, P, grad_P(i, :))
            ! formattage pour gnuplot
            plot_vec(i, :) = (/ x(i, :), (grad_P(i, :) / fnorme2(grad_P(i, :))) * 0.5_rp * dx /)
        end do
        call writeMat(plot_vec, "../sorties/grad_P.dat")


    case default
        write (*, *) "Mauvais choix de scénario."
    end select





    ! ---------------------------------------------------------------------------------------------------
    ! tension de surface Akinci / correction dimensionnelle (solver SPH)
    ! ---------------------------------------------------------------------------------------------------
    call normale(R, x, w, nor)
    do i = 1, np
        call F_TS(gamma_eau_air, i, x, w, nor, R, fts(i, :))
        nvec(i, :) = (/ x(i, :), (fts(i, :) / fnorme2(fts(i, :))) * 0.5_rp * dx /)
        !nvec(i, :) = (/ x(i, :), (nor(i, :) / fnorme2(nor(i, :))) * 0.5_rp * dx /)
    end do
    call writeMat(nvec, "../sorties/n_fts.dat")

    ! doit valoir zéro car cas à l'équilibre
    d_rwu_dt(:, 1) = w * grad_P(:, 1)
    d_rwu_dt(:, 2) = w * grad_P(:, 2)
    d_rwu_dt = d_rwu_dt + fts

    ! doit approcher (0, 0)
    print *, "sum_i [ w_i GR_p(P)_i + (F_TS)_i ] =", sum(d_rwu_dt(:, 1)), sum(d_rwu_dt(:, 2))

    ! champ de vecteur
    do i = 1, np
        plot_vec(i, :) = (/ x(i, :), d_rwu_dt(i, :) /)
    end do
    call writeMat(plot_vec, "../sorties/zero.dat")


    ! ---------------------------------------------------------------------------------------------------
    ! tension de surface précédente avec seulement le terme de cohésion
    ! ---------------------------------------------------------------------------------------------------
    !do i = 1, np
    !    call F_TS_cohesion(gamma_eau_air, i, x, w, R, fts(i, :))
    !    nvec(i, :) = (/ x(i, :), (fts(i, :) / fnorme2(fts(i, :))) * 0.5_rp * dx /)
    !    !nvec(i, :) = (/ x(i, :), (nor(i, :) / fnorme2(nor(i, :))) * 0.5_rp * dx /)
    !end do
    !d_rwu_dt(:, 1) = w * grad_P(:, 1)
    !d_rwu_dt(:, 2) = w * grad_P(:, 2)
    !d_rwu_dt = d_rwu_dt + fts
    !do i = 1, np
    !    plot_vec(i, :) = (/ x(i, :), d_rwu_dt(i, :) /)
    !end do
    !call writeMat(plot_vec, "../sorties/zero.dat")





    ! *******************************************************************************************************
    if (allocated(x1)) then
        deallocate(x1)
        deallocate(x2)
    end if
    deallocate(x)
    deallocate(w)
    deallocate(grad_P)
    deallocate(fts)
    deallocate(nor)
    deallocate(d_rwu_dt)

END PROGRAM main
