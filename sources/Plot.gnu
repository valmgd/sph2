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
plot "../sorties/x_points.dat"    u 1:2     lc rgb "#008000" lw 3 title "particules",\
     "../sorties/x_enveloppe.dat" u 1:2 w l lt rgb "green"   lw 1 title "bord de {/Symbol W}"




# ===========================================================================================================
# affichage écran
set term x11 enhanced
replot
