# Trois noyaux de cohésion
set term postscript eps enhanced color solid size 3.5,2.62 font 'Helvetica,12';
set output "../graphes/cohesion_akinci.ps";
set encoding utf8;
set title "Noyau de cohésion de N. Akinci 2D"
set xrange [0:1]
#set yrange [-1:2]
set xlabel "q = ||x_i - x|| / R"
set ylabel ""
set grid
set border 19
plot "../sorties/Cakinci.dat" u 1:2 w l lc rgb "blue" lw 3 title "C(x_i - x)"

reset;
set term postscript eps enhanced color solid size 3.5,2.62 font 'Helvetica,12';
set output "../graphes/cohesion_liu.ps";
set encoding utf8;
set title "Noyau de cohésion de M. Liu"
set xrange [0:1]
#set yrange [-1:2]
set xlabel "q = ||x_i - x|| / R"
set ylabel "C"
set grid
set border 19
plot "../sorties/Cliu.dat" u 1:2 w l lc rgb "blue" lw 3 title "C(x_i - x)", \
     "../sorties/AW1.dat" u 1:2 w l lc rgb "green" title "AW_1", \
     "../sorties/BW2.dat" u 1:2 w l lc rgb "red" title "BW_2"

set term postscript eps enhanced color solid size 3.5,2.62 font 'Helvetica,12';
set output "../graphes/cohesion_new.ps";
set encoding utf8;
set title "Noyau de cohésion"
set xrange [0:1]
#set yrange [-1:2]
set xlabel "q = ||x_i - x|| / R"
set ylabel ""
set grid
set border 19
plot "../sorties/Cnew.dat" u 1:2 w l lc rgb "blue" lw 3 title "C(x_i - x)"



reset;
set term postscript eps enhanced color solid size 3.5,2.62 font 'Helvetica,12';
set output "../graphes/gradR.ps";
set encoding utf8;
set title "Noyau de cohésion"
#set yrange [-1:2]
set xlabel "||x_i - x||"
set ylabel ""
set size ratio -1;
set grid
set border 19
plot "../sorties/P.dat" u 1:2:3 with points pointtype 31 pointsize 0.2 lc rgb "blue" title "particules",\
         "../sorties/x_enveloppe.dat" u 1:2 w l lt rgb "blue" lw 2 title "bord de {/Symbol W}",\
         "../sorties/nor.dat" u 1:2:3:4 with vectors head filled lt rgb "black" title "direction G^R(R)(x)";
#        "../sorties/gradR.dat" u 1:2:3:4 with vectors head filled lt rgb "black" title "direction G^R(R)(x)";

# ===========================================================================================================
# réinitialisation des paramètres
reset;
# fichier de sortie
#set term postscript eps enhanced color solid size 3.5,2.62 font 'Helvetica,12';
set term postscript eps enhanced color solid size 7,5.24 font 'Helvetica,12';
set output "../graphes/fluide.ps";
set encoding utf8;

d = system("cat ../sorties/x_d.dat");
#print "dimension : ", d;

# paramètres
set title "Gradient régularisé de la fonction 1";
set grid;
set palette rgbformulae 22,13,-31;
set cblabel "pression"
set key outside;
set border 16;

# duo de couleurs
set linetype 1 linecolor rgb "green"
set linetype 2 linecolor rgb "red"



if (d == 2) {
    ntics = 2;
    stats "../sorties/x_scale.dat" using 1 name "x" nooutput;
    stats "../sorties/x_scale.dat" using 2 name "y" nooutput;
    set xtics x_max/ntics;
    set ytics y_max/ntics;

#stats "../sorties/P.dat" using 3 name "P" nooutput
#set cbrange [0:P_max];
#set cbrange [0:0.15];
    set xlabel "x";
    set ylabel "y";
    set size ratio -1;

# tracé
#    plot "../sorties/P.dat" u 1:2:3 with points pointtype 5 pointsize 1 palette title "particules",\
#        "../sorties/x_enveloppe.dat" u 1:2 w l lt rgb "blue" lw 2 title "bord de {/Symbol W}",\
#        "../sorties/plot_vec.dat" u 1:2:3:4:($5>0?1:2) with vectors title "D{/Symbol rw}u / Dt" lw 0.5 linecolor variable,\
#"../sorties/normale.dat" u 1:2:3:4 with vectors head filled lt rgb "black" title "direction G^R(1)(x)",\

    plot "../sorties/P.dat" u 1:2:3 with points pointtype 31 pointsize 1 lc rgb "blue" title "particules",\
        "../sorties/x_enveloppe.dat" u 1:2 w l lt rgb "blue" lw 2 title "bord de {/Symbol W}",\
        "../sorties/plot_vec.dat" u 1:2:3:4:($5>0?1:2) with vectors title "D{/Symbol rw}u / Dt" lw 2 linecolor variable,\
        ;

} else {
    ntics = 2;
    stats "../sorties/x_scale.dat" using 1 name "x" nooutput;
    stats "../sorties/x_scale.dat" using 2 name "y" nooutput;
    stats "../sorties/x_scale.dat" using 3 name "z" nooutput;
    set xtics x_max/ntics;
    set ytics y_max/ntics;
    set ztics z_max/ntics;

#set cbrange [0:0.15];
    set view equal xyz;

    splot "../sorties/P.dat" u 1:2:3:4 with points pointtype 5 pointsize 1 palette title "particules",\
        "../sorties/x_enveloppe.dat" u 1:2:3 w l lt rgb "green" lw 2 title "bord de {/Symbol W}",\
        "../sorties/normale.dat" u 1:2:3:4:5:6 with vectors head filled lt rgb "black" title "direction G^R(1)(x)",\
        ;
}




# ===========================================================================================================
# affichage écran
set term x11 enhanced;
replot;
