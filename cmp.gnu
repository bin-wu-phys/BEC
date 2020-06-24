#!/usr/bin/gnuplot
set log
set xrange [1e-5:0.01]
set yrange [0.01:0.2]
set term x11 0
plot "f00400Qs05Nc3Nf3N2400dtm7V3.2.1.dat" u 1:2 index 302:342:2 w lp,\
"f00400Qs05Nc3Nf3N1600dtm7V3.2.1.dat" u 1:2 index 302:342:2 w lp,\
"f00400Qs05Nc3Nf3N3200dtm7V3.2.1.dat" u 1:2 index 302:342:2 w lp

reset
set log x
set xrange [1e-5:0.01]
set yrange [0.175:0.197]
set term x11 1 
plot "f00400Qs05Nc3Nf3N2400dtm7V3.2.1.dat" u 1:2 index 342:472:10 w lp,\
"f00400Qs05Nc3Nf3N1600dtm7V3.2.1.dat" u 1:2 index 342:472:10 w lp,\
"f00400Qs05Nc3Nf3N3200dtm7V3.2.1.dat" u 1:2 index 342:472:10 w lp,\
"f00400Qs05Nc3Nf3N4800dtm7V3.2.1.dat" u 1:2 index 342:472:10 w lp

reset
set log x
set xrange [1e-5:0.01]
set yrange [0.44:0.5]
set term x11 2 
plot "f00400Qs05Nc3Nf3N2400dtm7V3.2.1.dat" u 1:3 index 342:472:10 w lp,\
"f00400Qs05Nc3Nf3N1600dtm7V3.2.1.dat" u 1:3 index 342:472:10 w lp,\
"f00400Qs05Nc3Nf3N3200dtm7V3.2.1.dat" u 1:3 index 342:472:10 w lp,\
"f00400Qs05Nc3Nf3N4800dtm7V3.2.1.dat" u 1:3 index 472 w l lw 3

#set yrange [0.192:0.197]
#set term x11 2
#set yrange [0.192:0.197]
#set term x11 2
#plot "f00400Qs05Nc3Nf3N2400dtm7V3.2.1.dat" u 1:2 index 500:600:4 w lp,\
#"f00400Qs05Nc3Nf3N1600dtm7V3.2.1.dat" u 1:2 index 500:600:4 w lp,\
#"f00400Qs05Nc3Nf3N3200dtm7V3.2.1.dat" u 1:2 index 500:600:4 w lp,\
#"f00400Qs05Nc3Nf3N2400dtm7V3.2.1.dat" u 1:2 index 600 w lp

#set yrange [0.187:0.1955]
#set term x11 3
#plot "f00400Qs05Nc3Nf3N2400dtm7V3.2.1.dat" u 1:2 index 600:700:4 w lp,\
#"f00400Qs05Nc3Nf3N1600dtm7V3.2.1.dat" u 1:2 index 600:700:4 w lp,\
#"f00400Qs05Nc3Nf3N3200dtm7V3.2.1.dat" u 1:2 index 600:700:4 w lp,\
#"f00400Qs05Nc3Nf3N2400dtm7V3.2.1.dat" u 1:2 index 700 w lp

set term x11 4
reset
plot "f00400Qs05Nc3Nf3N1600dtm7V3.2.1mac.dat" u 1:6 w p,\
"f00400Qs05Nc3Nf3N3200dtm7V3.2.1mac.dat" u 1:6 w p,\
"f00400Qs05Nc3Nf3N4000dtm7V3.2.1mac.dat" u 1:6 w p
pause -1

