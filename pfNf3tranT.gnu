#!/usr/bin/gnuplot
reset
set term post  portrait color enh dl 3.5 size 9.6,7.2
set out "pfNf3tranT.eps"
set size ratio .8

Qs = 0.5

set key font "Helvetica,25"
set xtics font "Helvetica,28" #offset 1
set ytics font "Helvetica,28" #offset -1
set xlabel font "Helvetica,28"
set ylabel font "Helvetica,28"
set title font "Helvetica,28"

set style line 2 lw 12 lt 8 ps 1.5 pt 4 lc rgb "red"
set style line 3 lw 12 lt 2 pt 7 ps 1.1 lc rgb "#ff6600"
set style line 4 lw 12 lt 5 pt 6 ps 5 lc rgb "dark-green"
set style line 5 lw 6 lt 1 pt 7 ps 1.1 lc rgb "blue"
set style line 1 lw 6 lt 3 pt 6 ps 2 lc rgb "#800080"
set style line 6 lw 6 lt 1 pt 7 ps 1.1 lc rgb "black"

set key top left at 0.01/Qs, 0.002/Qs
set key spacing 4 samplen 6
set log

set xrange [2e-4:2/Qs]
set yrange [1e-4/Qs:.2/Qs]
#set format y "%1.0t{/Symbol \264}10^{%T}"
#set ytics 0.05/Qs,0.025/Qs,0.2/Qs

#set title "N_f = 3.0, f_0 = 0.1" offset -1
set xlabel "p/Q_s" offset 1
set ylabel "pf/Q_s" offset -5

feq(x)=x/((exp(8.008*Qs*x+8.008*0.0178521)-1))

plot "f00260Qs05Nc3Nf3N800dtm6V3.2.1.dat" u ($1/Qs):($2/Qs) index 190 w l ls 2 ti "{/Symbol t} &{Q_s} = &{0}{/Symbol t}_c&{0}",\
"f00260Qs05Nc3Nf3N800dtm6V3.2.1.dat" u ($1/Qs):($2/Qs) index 306 w l ls 3 ti "{/Symbol t} Q_s = &{0}1.5",\
"f00260Qs05Nc3Nf3N800dtm6V3.2.1.dat" u ($1/Qs):($2/Qs) index 506 w l ls 4 ti "{/Symbol t} Q_s = &{0}2.5",\
feq(x) w l ls 5 noti,\
"f00260Qs05Nc3Nf3N800dtm6V3.2.1.dat" u ($1/Qs):($2/Qs) index 2508 every 10 w p ls 1 ti "{/Symbol t} Q_s = 12.5",\
"f00260Qs05Nc3Nf3N800dtm6V3.2.1.dat" u ($1/Qs):($2/Qs) index 2507 every 10 w p ls 1 noti,\
"f00260Qs05Nc3Nf3N800dtm6V3.2.1.dat" u ($1/Qs):($2/Qs) index 2506 every 2 w p ls 1 noti

reset
set term x11

