# ===========================================================================================================
# réinitialisation des paramètres
reset
# fichier de sortie
set term postscript eps enhanced color solid size 3.5,2.62 font 'Helvetica,12'
set output "../graphes/pave.ps"
set encoding utf8

# paramètres
set title "Gradient régularisé de la fonction 1"
set grid
set xlabel "x"
set ylabel "y"
set size ratio -1 # repère orthonormé
set xrange[-0.4:1.4]
set yrange[-0.4:1.4]
set border 16

# tracé
plot "../sorties/x.dat" u 1:2 lc rgb "#008000" lw 3 title "particules",\
     "../sorties/nvec.dat" u 1:2:3:4 with vectors head filled lt rgb "red" title "G^R(1)(x)",\
     "../sorties/cc.dat" u 1:2 w l lw 1 lt rgb "green" title "bord de {/Symbol W}"
#    "../sorties/grad.dat" u 1:2:3:4 with vectors head filled lt rgb "green"



# ===========================================================================================================
reset
set term postscript eps enhanced color solid size 3.5,2.62 font 'Helvetica,12'
set output "../graphes/circle.ps"
set encoding utf8

# paramètres
set title "Discrétisation"
set grid
set xlabel "x"
set ylabel "y"
set size ratio -1 # repère orthonormé
#set xrange[-0.2:2.2]
#set yrange[-0.2:2.2]
#set border 16

# tracé
plot "../sorties/x_circle.dat" u 1:2 lc rgb "#008000" lw 1 title "particules",\
     "../sorties/cc.dat" u 1:2 w l lw 1 title "bord de {/Symbol W}"




# ===========================================================================================================
reset
set term postscript eps enhanced color solid size 7,5.24 font 'Helvetica,12'
set output "../graphes/pression.ps"
set encoding utf8

# paramètres
set title "Répartition de la pression\ndans une bulle d'eau à l'équilibre"
set grid
set key outside
set key right top
set xlabel "x"
set ylabel "y"
set size ratio -1 # repère orthonormé
#set xrange[-0.5:2.5]
#set yrange[-0.5:2.5]
set palette rgbformulae 22,13,-31
set cblabel "pression"

# tracé
plot "../sorties/P.dat" u 1:2:3 with points pointtype 5 pointsize 1 palette title "particule",\
     "../sorties/cc.dat" u 1:2 w l lw 1 title "délimitation bulle",\
     "../sorties/zero.dat" u 1:2:3:4 with vectors head filled lt rgb "black" title "wi * GR_p(P) + F_TS"
#     "../sorties/n_fts.dat" u 1:2:3:4 with vectors head filled lt rgb "green" title "Force tension surface"


#    "../sorties/grad_P.dat" u 1:2:3:4 with vectors head filled lt rgb "red" title "gradient de pression",\

# ===========================================================================================================
# réinitialisation des paramètres
reset
# fichier de sortie
set term postscript eps enhanced color solid size 3.5,2.62 font 'Helvetica,12'
set output "../graphes/Cakinci.ps"
set encoding utf8

# paramètres
set title "Force de cohésion de N. Akinci"
set grid
set xlabel "||x_i - x||"
set ylabel "C"
set border 3

# tracé
plot "../sorties/CAkinci.dat" u 1:2 w l lt rgb "blue" lw 1 title "C(||x_i - x||)"



# ===========================================================================================================
# affichage écran
set term x11 enhanced
replot
