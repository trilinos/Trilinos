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
set y2range [1.0e-13:1]
set logscale x
set logscale y2

plot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-2.dat" using 1:2 w l ls 1 title "Exact"
replot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-2.dat" using 1:3 w l ls 2 title "Backward Euler"
replot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-2.dat" using 1:5 w l ls 3 axes x1y2 title "dt"
#set title ""
set xlabel "Time"
set ylabel "Solution"
set ytics nomirror
set y2tics nomirror
set y2label "dt" offset -3,0
replot
set terminal postscript eps solid color lw 1
set output "BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-2.eps"
replot
set term x11

plot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-3.dat" using 1:2 w l ls 1 title "Exact"
replot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-3.dat" using 1:3 w l ls 2 title "Backward Euler"
replot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-3.dat" using 1:5 w l ls 3 axes x1y2 title "dt"
#set title ""
set xlabel "Time"
set ylabel "Solution"
set ytics nomirror
set y2tics nomirror
set y2label "dt" offset -3,0
replot
set terminal postscript eps solid color lw 1
set output "BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-3.eps"
replot
set term x11

plot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-4.dat" using 1:2 w l ls 1 title "Exact"
replot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-4.dat" using 1:3 w l ls 2 title "Backward Euler"
replot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-4.dat" using 1:5 w l ls 3 axes x1y2 title "dt"
#set title ""
set xlabel "Time"
set ylabel "Solution"
set ytics nomirror
set y2tics nomirror
set y2label "dt" offset -3,0
replot
set terminal postscript eps solid color lw 1
set output "BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-4.eps"
replot
set term x11

plot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-5.dat" using 1:2 w l ls 1 title "Exact"
replot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-5.dat" using 1:3 w l ls 2 title "Backward Euler"
replot "../../../../../build/packages/rythmos/test/ConvergenceTest/BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-5.dat" using 1:5 w l ls 3 axes x1y2 title "dt"
#set title ""
set xlabel "Time"
set ylabel "Solution"
set ytics nomirror
set y2tics nomirror
set y2label "dt" offset -3,0
replot
set terminal postscript eps solid color lw 1
set output "BackwardEuler_FirstOrderError_var_dt_RelError\=1.0e-5.eps"
replot
set term x11
