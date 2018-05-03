# ===========================================================================================================
# réinitialisation des paramètres
reset;
# fichier de sortie
set term postscript eps enhanced color solid size 3.5,2.62 font 'Helvetica,12';
set output "../graphes/pave.ps";
set encoding utf8;

d = system("cat ../sorties/x_d.dat");
print "dimension : ", d;

# paramètres
set title "Gradient régularisé de la fonction 1";
set grid;

if (d == 2) {
    set xlabel "x";
    set ylabel "y";
    set border 16;
    set size ratio -1;

# tracé
    plot "../sorties/x_points.dat" u 1:2 lc rgb "#008000" lw 2 title "particules",\
        "../sorties/x_enveloppe.dat" u 1:2 w l lt rgb "green" lw 2 title "bord de {/Symbol W}";
} else {
    set view equal xyz
    splot "../sorties/x_points.dat" u 1:2:3 lc rgb "#008000" lw 2 title "particules",\
        "../sorties/x_enveloppe.dat" u 1:2:3 w l lt rgb "green" lw 2 title "bord de {/Symbol W}";
}




# ===========================================================================================================
# affichage écran
set term x11 enhanced;
replot;
