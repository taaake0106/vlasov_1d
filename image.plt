reset 
set term png
set output "epot.jpg"
set xrange [10:1000]
set yrange [1e-8:1]
set title "e pot"
set xlabel "N_m"
set ylabel "relative error"
set size square
set logscale x 
set logscale y 
set format x "10^{%L}"
set format y "10^{%L}"
set grid
set key left bottom
plot "../read/e_pot_2nd.dat" using 1:2 with lp title "2nd", "../read/e_pot_4th.dat" using 1:2 with lp title "4th", "../read/e_pot_6th.dat" using 1:2 with lp title "6th",100/(x*x) with lines,10000/(x*x*x*x) with lines,1000000/(x*x*x*x*x*x) with lines