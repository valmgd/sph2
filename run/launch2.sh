#!/bin/bash

rayon='0.0001 0.0010 0.0100 0.1000 1.0000'
sigma='0.05 0.06 0.07 0.08 0.09'
intervalle='2 3 4 5 6'

rm -f rayon.dat
rm -f sigma.dat
rm -f intervalle.dat

for elt in $rayon; do
    echo '#dimension' > ../entrees/constantes
    echo '2' >> ../entrees/constantes
    echo '#sigma' >> ../entrees/constantes
    echo '0.07D0' >> ../entrees/constantes
    echo '#intervalle' >> ../entrees/constantes
    echo '4.0D0' >> ../entrees/constantes
    echo '#n' >> ../entrees/constantes
    echo '11' >> ../entrees/constantes
    echo '#bornes' >> ../entrees/constantes
    echo '0.0D0 1.0D0' >> ../entrees/constantes
    echo '0.0D0 1.0D0' >> ../entrees/constantes
    echo '0.0D0 1.0D0' >> ../entrees/constantes
    echo '#centre' >> ../entrees/constantes
    echo '0.0D0 0.0D0 0.0D0' >> ../entrees/constantes
    echo '#rayon' >> ../entrees/constantes
    echo "$elt" >> ../entrees/constantes

    ./truc
    (cd ../sources; gnuplot Plot.gnu)
    (cd ../graphes; mv fluide.ps "fluide_rayon$elt.eps")
    echo "###"
    cat ../entrees/constantes
    echo
done

for elt in $sigma; do
    echo '#dimension' > ../entrees/constantes
    echo '2' >> ../entrees/constantes
    echo '#sigma' >> ../entrees/constantes
    echo "$elt" >> ../entrees/constantes
    echo '#intervalle' >> ../entrees/constantes
    echo '4.0D0' >> ../entrees/constantes
    echo '#n' >> ../entrees/constantes
    echo '11' >> ../entrees/constantes
    echo '#bornes' >> ../entrees/constantes
    echo '0.0D0 1.0D0' >> ../entrees/constantes
    echo '0.0D0 1.0D0' >> ../entrees/constantes
    echo '0.0D0 1.0D0' >> ../entrees/constantes
    echo '#centre' >> ../entrees/constantes
    echo '0.0D0 0.0D0 0.0D0' >> ../entrees/constantes
    echo '#rayon' >> ../entrees/constantes
    echo '0.01' >> ../entrees/constantes

    ./truc
    (cd ../sources; gnuplot Plot.gnu)
    (cd ../graphes; mv fluide.ps "fluide_sigma$elt.eps")
    echo "###"
    cat ../entrees/constantes
    echo
done

for elt in $intervalle; do
    echo '#dimension' > ../entrees/constantes
    echo '2' >> ../entrees/constantes
    echo '#sigma' >> ../entrees/constantes
    echo '0.07D0' >> ../entrees/constantes
    echo '#intervalle' >> ../entrees/constantes
    echo "$elt" >> ../entrees/constantes
    echo '#n' >> ../entrees/constantes
    echo '11' >> ../entrees/constantes
    echo '#bornes' >> ../entrees/constantes
    echo '0.0D0 1.0D0' >> ../entrees/constantes
    echo '0.0D0 1.0D0' >> ../entrees/constantes
    echo '0.0D0 1.0D0' >> ../entrees/constantes
    echo '#centre' >> ../entrees/constantes
    echo '0.0D0 0.0D0 0.0D0' >> ../entrees/constantes
    echo '#rayon' >> ../entrees/constantes
    echo '0.01' >> ../entrees/constantes

    ./truc
    (cd ../sources; gnuplot Plot.gnu)
    (cd ../graphes; mv fluide.ps "fluide_intervalle$elt.eps")
    echo "###"
    cat ../entrees/constantes
    echo
done

cp ../entrees/backup ../entrees/constantes
