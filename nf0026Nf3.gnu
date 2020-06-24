#!/usr/bin/gnuplot
reset
set term post portrait color enh dl 3.5 size 8,6
set out "nf0026Nf3.eps"
set size ratio .8

Qs = 0.5
Qs3 = Qs*Qs*Qs

set key spacing 3.5 samplen 6
#set format y "%1.0t{/Symbol \264}10^{%T}"
set ytics 0.070234,0.000004,0.070250

set key font "Helvetica,22"
set tics font "Helvetica,25"
set xlabel font "Helvetica,25"
set ylabel font "Helvetica,25"
set title font "Helvetica,25"

set key center center

set style line 1 lw 8 lt 1 pt 7 ps 1.1 lc rgb "red"
set style line 2 lw 12 lt 1 pt 7 ps 1.1 lc rgb "#ff6600"
set style line 3 lw 12 lt 2 pt 7 ps 1.1 lc rgb "#d38d5f"
set style line 4 lw 12 lt 5 pt 7 ps 1.1 lc rgb "dark-green"
set style line 5 lw 12 lt 3 pt 7 ps 2 lc rgb "red"

set xrange [0:4.0*Qs]
set yrange [0.070234:0.070250]
#set format y "%1.0t{/Symbol \264}10^{%T}"

set xlabel "{/Symbol t} Q_s" offset 1
set ylabel "n/Q_s^3" offset -3

plot "f00260Qs05Nc3Nf3N800dtm6V3.2.1mac.dat" u ($1*Qs):($6/Qs3) w l ls 1 noti

reset
set term x11

