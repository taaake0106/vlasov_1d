reset
set term gif animate
set output "MPP7_kk_j=0.5.gif"
set xrange[0.0:1.0]
set yrange[-1.0:1.0]
set pm3d 
set pm3d map
set ticslevel 0
set cbrange[0:3.1]
set xlabel "x/L"
set ylabel "v/V"
set size square
set palette rgbformulae 33,13,10

do for [i=0:2560]{
    if(i%10==0){
set title sprintf("t=%d",i)
     splot  sprintf("../SL-MPP7_%d.dat",i) with pm3d

    }
}

