touch plotcommand.gnu
TIMETITLE="Runtime vs. Number Of Nodes for Matrix-Matrix Multiply (lower is better)"
TIMEOUT="time.ps"
EFFTITLE="Scale efficiency vs. Number Of Nodes for Matrix-Matrix Multiply (higher is better)"
EFFOUT="eff.ps" 
echo "set term postscript color  \\" >> plotcommand.gnu
echo "dashed dashlength 3.0 linewidth 3.0 " >> plotcommand.gnu
echo "set yrange [0:]" >> plotcommand.gnu
echo "set data style linespoints" >> plotcommand.gnu
echo "set style line 2 linecolor rgb \"orange\" " >> plotcommand.gnu
echo "set style line 3 linecolor rgb \"red\" " >> plotcommand.gnu
echo "set style line 4 linecolor rgb \"purple\" " >> plotcommand.gnu
echo "set xlabel 'Number of Nodes'" >> plotcommand.gnu
echo "set log y" >> plotcommand.gnu
echo "set yrange [.1:]" >> plotcommand.gnu
echo "set key center right" >> plotcommand.gnu
if [ $1 == "time" ]
then
echo "set output '$TIMEOUT'" >> plotcommand.gnu
echo "set title '$TIMETITLE'" >> plotcommand.gnu
echo "set ylabel 'Runtime (seconds)'" >> plotcommand.gnu
echo "plot 'ttimes.out' ls 4 title 'Tpetra Matrix-Matrix Multiply (with all optimizations)', 'etimes.out' ls 3 title 'EpetraExt Matrix-Matrix Multiply'" >> plotcommand.gnu
else
echo "set title '$EFFTITLE'" >> plotcommand.gnu
echo "set output '$EFFOUT'" >> plotcommand.gnu
echo "set ylabel 'Percent Efficiency'" >> plotcommand.gnu
echo "plot 'teffs.out' ls 4 title 'Tpetra Matrix-Matrix Multiply', 'eeffs.out' ls 3 title 'EpetraExt Matrix-Matrix Multiply'" >> plotcommand.gnu
fi
gnuplot plotcommand.gnu
rm -f plotcommand.gnu
if [ $1 == "time" ]
then
convert -density 150 -geometry 100% -rotate 90 time.ps time.jpg
else
convert -density 150 -geometry 100% -rotate 90 eff.ps eff.jpg
fi
rm *.ps
