#!/usr/bin/gnuplot
reset
set xrange [1e-5:0.1]
set log
plot "f00400Qs05Nc3Nf3N4800dtm8V3.2.1.dat" u 1:(abs((0.177749676603747-0.5*$1)/$2-1)) index 342 ti "pf",\
"f00400Qs05Nc3Nf3N800dtm6V3.2.1.dat" u 1:(abs((0.177749676603747-0.5*$1)/$2-1)) index 342 ti "pf",\
"f00400Qs05Nc3Nf3N4800dtm8V3.2.1.dat" u 1:(abs((0.5-1.42259*$1)/$3-1))  index 342 ti "F",\
"f00400Qs05Nc3Nf3N800dtm6V3.2.1.dat" u 1:(abs((0.5-1.42259*$1)/$3-1))  index 342 ti "F",\
"f00400Qs05Nc3Nf3N4800dtm8V3.2.1.dat" u 1:(abs((0.000109155+ 1.08036*$1*$1)/$4-1))  index 342 ti "F_f",\
"f00400Qs05Nc3Nf3N800dtm6V3.2.1.dat" u 1:(abs((0.000109155+ 1.08036*$1*$1)/$4-1))  index 342 ti "F_f"
pause -1
