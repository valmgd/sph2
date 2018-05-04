#!/bin/bash

rayon='0.0001 0.0010 0.0100 0.1000 1.0000'
sigma='0.5 0.6 0.7 0.8 0.9'
intervalle='2 3 4 5 6'

rm -f rayon.dat
rm -f sigma.dat
rm -f intervalle.dat

for elt in $rayon; do
    echo '#dimension' > ../entrees/constantes
    echo '2' >> ../entrees/constantes
    echo '#sigma' >> ../entrees/constantes
    echo '0.073D0' >> ../entrees/constantes
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

    texte=`./truc`
    echo "$elt     $texte" >> rayon.dat
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
    echo '0.7' >> ../entrees/constantes

    texte=`./truc`
    echo "$elt     $texte" >> sigma.dat
done

for elt in $intervalle; do
    echo '#dimension' > ../entrees/constantes
    echo '2' >> ../entrees/constantes
    echo '#sigma' >> ../entrees/constantes
    echo '0.073D0' >> ../entrees/constantes
    echo '#intervalle' >> ../entrees/constantes
    echo "$intervalle" >> ../entrees/constantes
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

    texte=`./truc`
    echo "$elt     $texte" >> intervalle.dat
done

cp ../entrees/backup ../entrees/constantes
