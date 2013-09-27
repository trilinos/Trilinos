set size 0.6,0.6
set style line  1 linetype  1  linewidth 2
set style line  2 linetype  2  linewidth 2
set style line  3 linetype  3  linewidth 2
set style line  4 linetype  4  linewidth 2
set style line  5 linetype  5  linewidth 2
set style line  6 linetype  6  linewidth 2
set format x "%g"
set format y "%g"

set xrange [1.0e-12:1]
set yrange [0.0:1]
set logscale x
plot "logtime.dat" using 1:2 w l ls 1 title "x"
set title "Log-Time Problem"
set xlabel "Time"
set ylabel "Solution"
replot
set terminal postscript eps solid color lw 1
set output "logtime-log.eps"
replot
set term x11

set xrange [1.0e-12:1]
set yrange [0.0:1]
set nologscale x
plot "logtime.dat" using 1:2 w l ls 1 title "x"
set title "Log-Time Problem"
set xlabel "Time"
set ylabel "Solution"
replot
set terminal postscript eps solid color lw 1
set output "logtime-linear.eps"
replot
set term x11
