set size 0.6,0.6
set style line  1 linetype  1  linewidth 2
set style line  2 linetype  2  linewidth 2
set style line  3 linetype  3  linewidth 2
set style line  4 linetype  4  linewidth 2
set style line  5 linetype  5  linewidth 2
set style line  6 linetype  6  linewidth 2
set format y "%g"

#set xrange [0.2:0.5]
plot "sincos.dat" using 1:2 w l ls 1 title "x0"
replot "sincos.dat" using 1:3 w l ls 2 title "x1"
set title "Sine-Cosine Problem"
set xlabel "Time"
set ylabel "Solution"
replot
set terminal postscript eps solid color lw 1
set output "sincos.eps"
replot
set term x11
