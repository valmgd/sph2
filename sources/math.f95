! ===========================================================================================================
! Module contenant les fonctions :
!
! linspace
! seq
! evaluate
! evaluate2D
! saveSol
! saveMat
! affMat
! id
! plu
! cholesky
! prodScal
! prodMatVec
! descenteRemonte
! linSolve
! norme2
! fnorme2
! normeInf
! trapeze
! trieBulle
! readInteger
! readReal
! newton1D
! newtonND
! ===========================================================================================================

MODULE math

    implicit none

    integer, parameter :: sp = kind(1.0E0)
    integer, parameter :: dp = kind(1.0D0)
    integer, parameter :: rp = dp
    real(rp), parameter :: pi = 4.0_rp * atan(1.0_rp)

contains

    ! -------------------------------------------------------------------------------------------------------
    ! créer un vecteur uniforme en précisant la taille
    ! -------------------------------------------------------------------------------------------------------
    ! a : valeur de début
    ! b : valeur de fin
    ! n : nombre de valeur du vecteur retourné
    function linspace(a, b, n)
        ! param
        real(rp), intent(in) :: a, b
        integer, intent(in) :: n

        ! return
        real(rp), dimension(n) :: linspace

        ! variables locales
        real(rp) :: h
        integer :: i

        h = (b - a) / (n - 1)
        !linspace(1:n) = (a + i * h, i = 0, n - 1)
        do i = 0, n - 2
            linspace(i + 1) = a + i * h
        end do
        linspace(n) = b
    end function linspace



    ! -------------------------------------------------------------------------------------------------------
    ! créer un vecteur uniforme en précisant le pas
    ! -------------------------------------------------------------------------------------------------------
    ! a : valeur de début
    ! b : valeur de fin
    ! h : pas éventuellement diminué pour que la valeur finale soit b
    function seq(a, b, h)
        ! param
        real(rp), intent(in) :: a, b
        real(rp) :: h

        ! return
        real(rp), dimension(:), allocatable :: seq

        ! variables locales
        integer :: n, i

        ! nombre de valeurs
        n = ceiling((b - a) / h) + 1
        allocate(seq(n))

        ! nouveau pas inf ou égal au param
        h = (b - a) / (n - 1)
        seq = linspace(a, b, n)
    end function seq



    ! -------------------------------------------------------------------------------------------------------
    ! evalue une fonction sur un vecteur
    ! -------------------------------------------------------------------------------------------------------
    ! fc : fonction de une variable
    ! vec : vecteur quelconque
    ! n : taille de vec
    !
    ! sortie : fc(vec) qui est un vecteur de taille n
    function evaluate(fc, vec)
        ! param
        real(rp), external :: fc
        real(rp), dimension(:), intent(in) :: vec

        ! return
        real(rp), dimension(size(vec)) :: evaluate

        ! variables locales
        integer :: i

        do i = 1, size(vec)
            evaluate(i) = fc(vec(i))
        end do
    end function

    function evaluate2D(fc, vx, vy)
        ! PARAMÈTRES
        real(rp), external :: fc
        real(rp), dimension(:), intent(in) :: vx
        real(rp), dimension(:), intent(in) :: vy

        ! RETURN
        real(rp), dimension(size(vx)) :: evaluate2D

        ! VARIABLES LOCALES
        integer :: i

        do i = 1, size(vx)
            evaluate2D(i) = fc(vx(i), vy(i))
        end do
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! sauvegarde les antécédants et la solution en colonne dans un fichier
    ! -------------------------------------------------------------------------------------------------------
    ! x : abscisse
    ! sol : ordonnée
    ! n : taille de x et sol
    ! fichier : nom du fichier créé
    subroutine saveSol(x, sol, fichier)
        ! parametres
        real(rp), dimension(:) :: x, sol
        character(len=*) :: fichier

        ! variables locales
        integer :: i

        open(unit=0, file=fichier)
        do i = 1, size(x)
            write (0, *) x(i), sol(i)
        end do
        close(0)
    end subroutine saveSol



    ! -------------------------------------------------------------------------------------------------------
    ! sauvegarde les antécédants et la solution en colonne dans un fichier
    ! -------------------------------------------------------------------------------------------------------
    ! x : abscisse
    ! mat : matrice de n lignes
    ! n : taille de x
    ! fichier : nom du fichier créé
    subroutine saveMat(x, mat, n, fichier)
        ! parametres
        real(rp), dimension(:) :: x
        real(rp), dimension(:, :) :: mat
        integer :: n
        character(len=*) :: fichier

        ! variables locales
        integer :: i

        open(unit=0, file=fichier)
        do i = 1, n
            write (0, *) x(i), mat(i, :)
        end do
        close(0)
    end subroutine saveMat



    ! -------------------------------------------------------------------------------------------------------
    ! affichage matrice en console
    ! -------------------------------------------------------------------------------------------------------
    ! n : nombre de ligne de la matrice
    ! mat : matrice à afficher
    subroutine affMat(mat)
        ! parametres
        real(rp), dimension(:, :), intent(in) :: mat

        ! variables locales
        integer :: i, n
        integer, dimension(2) :: forme

        forme = shape(mat)
        n = forme(1)

        do i = 1, n
            print '(100F7.1)', mat(i, :)
        end do
    end subroutine affMat



    ! -------------------------------------------------------------------------------------------------------
    ! matrice identité
    ! -------------------------------------------------------------------------------------------------------
    ! n : taille de la matrice à créer
    ! mat : matrice identité n par n retournée
    subroutine id(n, mat)
        ! parametres
        integer, intent(in) :: n
        real(rp), dimension(n, n), intent(out) :: mat

        ! variables locales
        integer :: i

        mat = 0.d0
        do i = 1, n
            mat(i, i) = 1.d0
        end do
    end subroutine id



    ! -------------------------------------------------------------------------------------------------------
    ! décomposition PLU
    ! -------------------------------------------------------------------------------------------------------
    ! n : dimension de l'espace
    ! A : matrice à décomposer
    !
    ! retours :
    ! P : matrice d'interversion ligne si pivot nul
    ! L : triang inf
    ! U : triang sup
    subroutine plu(n, A, P, L, U)
        ! parametres
        integer, intent(in) :: n
        real(rp), dimension(n, n), intent(in) :: A
        real(rp), dimension(n, n), intent(out) :: P, L, U

        ! variables locales
        integer :: i, j, k
        real(rp), dimension(n) :: temp
        real(rp) :: maxk

        call id(n, P)
        call id(n, L)
        U = A

        do k = 1, n - 1

            ! détermination de i tq A(i, k) = max_j(A(j, k)) j>=k
            i = k
            maxk = abs(A(i, k))
            do j = k + 1, n
                if (abs(A(j, k)) > maxk) then
                    maxk = abs(A(j, k))
                    i = j
                end if
            end do

            ! permutation lignes colonnes
            if (i /= k) then
                ! échange des lignes k et i de P, L, U
                P((/k, i/), :) = P((/i, k/), :)
                L((/k, i/), :) = L((/i, k/), :)
                U((/k, i/), :) = U((/i, k/), :)

                ! échange des colonnes k et i de L
                L(:, (/k, i/)) = L(:, (/i, k/))
            end if

            do i = k + 1, n
                L(i, k) = U(i, k) / U(k, k)
                do j = k + 1, n
                    U(i, j) = U(i, j) - L(i, k) * U(k, j)
                end do
                U(i, k) = 0.d0
            end do

        end do
    end subroutine plu



    ! -------------------------------------------------------------------------------------------------------
    ! décomposition Cholesky (cas A SDP)
    ! -------------------------------------------------------------------------------------------------------
    ! n : dimension d'espace
    ! A : matrice à décomposer
    ! L : matrice symétrique issue de la décomposition (n * n)
    subroutine cholesky(n, A, L)
        ! paramètres
        integer, intent(in) :: n
        real(rp), dimension(n, n), intent(in) :: A
        real(rp), dimension(n, n), intent(out) :: L

        ! variables locales
        integer :: i, j

        L(1, 1) = sqrt(A(1, 1))

        do i = 2, n
            L(i, 1) = A(i, 1) / L(1, 1)
            L(i, i) = sqrt(A(i, i) - sum(L(i, 1:i-1)**2))
            do j = i + 1, n
                L(j, i) = (A(j, i) - sum(L(i, 1:i-1) * L(j, 1:i-1))) / L(i, i)
            end do
        end do
    end subroutine cholesky



    ! -------------------------------------------------------------------------------------------------------
    ! produit vecteur vecteur
    ! -------------------------------------------------------------------------------------------------------
    ! n : dimension de l'espace
    ! vec1, vec2 : deux vecteurs de taille n
    ! prod : scalaire, résultat du produit
    subroutine prodScal(n, vec1, vec2, prod)
        ! parametres
        integer, intent(in) :: n
        real(rp), dimension(n), intent(in) :: vec1, vec2
        real(rp), intent(out) :: prod

        ! variables locales
        integer :: i

        prod = 0.d0

        do i = 1, n
            prod = prod + vec1(i) * vec2(i) 
        end do
    end subroutine prodScal



    ! -------------------------------------------------------------------------------------------------------
    ! produit matrice vecteur
    ! -------------------------------------------------------------------------------------------------------
    ! n : dimension de l'espace
    ! mat : matrice n*n
    ! vec : vecteur de taille n
    ! prod : vecteur de taille n résultat du produit
    subroutine prodMatVec(n, mat, vec, prod)
        ! parametres
        integer, intent(in) :: n
        real(rp), dimension(n, n), intent(in) :: mat
        real(rp), dimension(n), intent(in) :: vec
        real(rp), dimension(n), intent(out) :: prod

        ! variables locales
        integer :: i
        real(rp) :: temp

        do i = 1, n
            call prodScal(n, mat(i, :), vec, temp)
            prod(i) = temp
        end do
    end subroutine prodMatVec



    ! -------------------------------------------------------------------------------------------------------
    ! résolution x tq Ly = Pb avec Ux = y
    ! -------------------------------------------------------------------------------------------------------
    ! n : dimension de l'espace
    ! b : second membre
    ! P : passage
    ! L : triang ing
    ! U : triang sup
    ! x : solution retournée
    subroutine descenteRemontee(n, b, P, L, U, x)
        ! parametres
        integer, intent(in) :: n! dimension de l'espace
        real(rp), dimension(n), intent(in) :: b! second menbre
        real(rp), dimension(n, n), intent(in) :: P, L, U! issu décomposition plu
        real(rp), dimension(n), intent(out) :: x! solution

        ! variables lodales
        integer :: i
        real(rp), dimension(n) :: Pb, verif
        real(rp) :: temp
        real(rp), dimension(n) :: y

        call prodMatVec(n, P, b, Pb)

        ! algorithme de desente y tq Ly = Pb
        y(1) = Pb(1) / L(1, 1)

        do i = 2, n
            call prodScal(i - 1, y(1:i-1), L(i, 1:i-1), temp)
            y(i) = (Pb(i) - temp) / L(i, i)
        end do

        ! algorithme de remontée x tq Ux = y
        !print *, U(n, n), y(n)
        x(n) = y(n) / U(n, n)

        do i = n - 1, 1, -1
            call prodScal(n - i, x(i+1:n), U(i, i+1:n), temp)
            x(i) = (y(i) - temp) / U(i, i)
        end do

        call prodMatVec(n, L, y, verif)
    end subroutine descenteRemontee



    ! -------------------------------------------------------------------------------------------------------
    ! résolution Ax = b
    ! -------------------------------------------------------------------------------------------------------
    ! method : au choix parmi plu, cholesky
    ! n : dimension de l'espace
    ! A : matrice n*n
    ! b : second membre
    ! x : solution retournée
    subroutine linSolve(method, n, A, b, x)
        ! parametres
        character(len = *), intent(in) :: method! plu ou cholesky
        integer, intent(in) :: n! dimension de l'espace
        real(rp), dimension(n, n), intent(in) :: A! application linéaire
        real(rp), dimension(n), intent(in) :: b! second membre
        real(rp), dimension(n), intent(out) :: x! solution

        ! variables locales
        real(rp), dimension(n, n) :: P, L, U

        if (method == "plu") then
            call plu(n, A, P, L, U)
        else if (method == "cholesky") then
            call id(n, P)
            call cholesky(n, A, L)
            U = transpose(L)
        else
            print *, "choix de méthode invalide"
        end if

        call descenteRemontee(n, b, P, L, U, x)
    end subroutine linSolve



    ! -------------------------------------------------------------------------------------------------------
    ! norme euclidienne discrète
    ! -------------------------------------------------------------------------------------------------------
    ! v : vecteur
    ! norme : valeur retournée
    subroutine norme2(v, norme)
        ! paramètres
        real(rp), dimension(:), intent(in) :: v
        real(rp), intent(out) :: norme

        norme = sqrt(sum(v**2))
    end subroutine norme2

    function fnorme2(v) result(norme)
        ! paramètres
        real(rp), dimension(:), intent(in) :: v

        ! return
        real(rp) :: norme

        norme = sqrt(sum(v**2))
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! norme infinie discrète
    ! -------------------------------------------------------------------------------------------------------
    ! v : vecteur
    ! norme : valeur retournée
    subroutine normeInf(v, norme)
        ! paramètres
        real(rp), dimension(:), intent(in) :: v
        real(rp), intent(out) :: norme

        norme = maxval(abs(v))
    end subroutine normeInf



    ! -------------------------------------------------------------------------------------------------------
    ! méthode des trapèzes
    ! -------------------------------------------------------------------------------------------------------
    ! x : vecteur des antécédants (subdivision abscisse)
    ! fx : vecteur des images (fx = f(x))
    ! somme : approximation de l'intégrale retournée
    subroutine trapeze(x, fx, somme)
        ! paramètres
        real(rp), dimension(:), intent(in) :: x, fx
        real(rp), intent(out) :: somme

        ! variables locales
        integer :: i, n

        somme = 0.d0
        n = size(x)

        do i = 1, n - 1
            somme = somme + (x(i + 1) - x(i)) * (fx(i) + fx(i + 1)) / 2
        end do
    end subroutine trapeze



    ! -------------------------------------------------------------------------------------------------------
    ! Trier un tableau avec l'aogo trie à bulle (peu performant)
    ! -------------------------------------------------------------------------------------------------------
    ! tab : tableau à trier
    subroutine trieBulle(tab)
        ! paramètres
        integer, dimension(:) :: tab

        ! variables locales
        integer :: i, j, temp, imax
        integer :: n

        n = size(tab)

        do j = n, 1, -1
            imax = 1
            temp = tab(1)
            do i = 1, j
                ! si on trouve un elt plus grand
                if (tab(i) > temp) then
                    imax = i
                    ! c'est le nouveau plus grand elt
                    temp = tab(i)
                end if
            end do

            ! on inverse le plus grand elt et le dernier elt
            temp = tab(imax)
            tab(imax) = tab(j)
            tab(j) = temp
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Parcours fichier à la recherche de balise de 2 caractères, lis nombre ligne suivante
    ! -------------------------------------------------------------------------------------------------------
    ! fichier : nom du fichier à lire
    ! balise : chaîne de deux caractères (ex : #n)
    function readInteger(fichier, balise)
        ! paramètres
        character(len = *), intent(in) :: fichier
        character(len = 2), intent(in) :: balise

        ! return
        integer :: readInteger

        ! variables locales
        character(len = 2) :: test = "--"

        open(unit = 1, file = fichier)
        read(1, *) test

        do while (test /= balise)
            read(1, *) test
        end do

        read(1, *) readInteger
        close(1)
    end function

    ! fichier : nom du fichier à lire
    ! balise : chaîne de deux caractères (ex : #n)
    function readReal(fichier, balise)
        ! paramètres
        character(len = *), intent(in) :: fichier
        character(len = 2), intent(in) :: balise

        ! return
        real(rp) :: readReal

        ! variables locales
        character(len = 2) :: test = "--"

        open(unit = 1, file = fichier)
        read(1, *) test

        do while (test /= balise)
            read(1, *) test
        end do

        read(1, *) readReal
        close(1)
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Méthode de Newton pour fonctions 1D (résolution f(x) = 0)
    ! -------------------------------------------------------------------------------------------------------
    ! x0 : itéré initial
    ! f : fonction dont on cherche le zéro
    ! fp : dérivée de f
    ! err : distance entre deux itérés en deçà de laquelle le code s'arrête
    ! itmax : si la distance n'est pas atteinte, nombre d'itération max
    ! x2 : sortie, approximation du zéro de f
    subroutine newton1D(x0, f, fp, err, itmax, x2)
        ! paramètres
        real(rp), intent(in) :: x0
        real(rp), external :: f, fp
        real(rp), intent(in) :: err
        integer, intent(in) :: itmax
        real(rp), intent(out) :: x2

        ! variables locales
        integer :: i
        real(rp) :: x1, dist

        ! pour être sûr de rentrer une fois dans le while
        dist = 1000000.0_rp
        x1 = x0
        i = 0

        do while ((dist > err) .and. (i <= itmax))
            x2 = x1 - f(x1) / fp(x1)
            dist = abs(x2 - x1)
            x1 = x2
            i = i + 1
        end do

        if (i == itmax) then
            write (*, *) "Program newton1D abort. Max iterations reach."
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Méthode de Newton pour fonctions ND (résolution f(x) = 0)
    ! -------------------------------------------------------------------------------------------------------
    ! x0 : itéré initial
    ! f : fonction dont on cherche le zéro (de R^n dans R^n) (voir interface)
    ! Jp : jacobienne de f (voir interface)
    ! err : distance entre deux itérés en deçà de laquelle le code s'arrête
    ! itmax : si la distance n'est pas atteinte, nombre d'itération max
    ! x2 : sortie, approximation du zéro de f
    subroutine newtonND(x0, f, Jf, err, itmax, x2)
        ! paramètres
        real(rp), dimension(:), intent(in) :: x0

        ! fonction de Rn dans Rn
        interface
            subroutine f(psi, s)
                real(8), dimension(:), intent(in) :: psi
                real(8), dimension(:), intent(out) :: s
            end subroutine
        end interface

        ! jacobienne de f
        interface
            subroutine Jf(psi, jacobienne)
                real(8), dimension(:), intent(in) :: psi
                real(8), dimension(:, :), intent(out) :: jacobienne
            end subroutine
        end interface

        real(rp), intent(in) :: err
        integer, intent(in) :: itmax
        real(rp), dimension(size(x0)), intent(out) :: x2

        ! variables locales
        real(rp), dimension(size(x0), size(x0)) :: A
        real(rp), dimension(size(x0)) :: b
        real(rp), dimension(size(x0)) :: temp
        integer :: i
        real(rp), dimension(size(x0)) :: x1
        real(rp) :: dist

        ! pour être sûr de rentrer une fois dans le while
        dist = err + 1.0_rp
        x1 = x0
        i = 0

        do while ((dist > err) .and. (i <= itmax))
            call Jf(x1, A)
            call f(x1, b)
            call linSolve("plu", size(x0), A, b, temp)
            x2 = x1 - temp
            call norme2(x2 - x1, dist)
            x1 = x2
            i = i + 1
        end do

        if (i == itmax) then
            write (*, *) "Program newton1D abort. Max iterations reach."
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Sauvegarde une matrice dans un fichier
    ! -------------------------------------------------------------------------------------------------------
    subroutine writeMat(mat, fichier)
        ! paramètres
        real(rp), dimension(:, :), intent(in) :: mat
        character(len=*), intent(in) :: fichier

        ! variables locales
        integer :: i
        integer, dimension(2) :: sh

        sh = shape(mat)

        open(unit = 10, file = fichier)
        do i = 1, sh(1)
            write (10, *) mat(i, :)
        end do
        close(10)
    end subroutine

END MODULE math
