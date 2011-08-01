touch plotcommand.gnu
TIMETITLE="Runtime vs. Number Of Nodes for Tpetra Matrix-Matrix Multiply (lower is better)"
TIMEOUT="time.ps"
EFFTITLE="Scale efficiency vs. Number Of Nodes for Tpetra Matrix-Matrix Multiply (higher is better)"
EFFOUT="eff.ps" 
echo "set term postscript color  \\" >> plotcommand.gnu
echo "dashed dashlength 3.0 linewidth 3.0 " >> plotcommand.gnu
echo "set key bottom right" >> plotcommand.gnu
echo "set style data linespoints" >> plotcommand.gnu
echo "set style line 4 linecolor rgb \"purple\" " >> plotcommand.gnu
echo "set style line 3 linecolor rgb \"red\" " >> plotcommand.gnu
echo "set style line 2 linecolor rgb \"blue\" " >> plotcommand.gnu
echo "set yrange [0:]" >> plotcommand.gnu
echo "set xlabel 'Number of Nodes'" >> plotcommand.gnu
if [ $1 == "time" ]
then
echo "set ylabel 'Runtime (seconds)'" >> plotcommand.gnu
echo "set title '$TIMETITLE'" >> plotcommand.gnu
echo "set output '$TIMEOUT'" >> plotcommand.gnu
echo "plot 'oldtimes.out' ls 3 title 'Initial Matrix-Matrix Multiply', 'newtimes.out' ls 2 title 'Matrix-Matrix Multiply With Optimzations From 3.1 and 3.2', 'new2times.out' ls 4 title 'Matrix-Matrix Multiply With All Optimizations" >> plotcommand.gnu
else
echo "set ylabel 'Percent Efficiency'" >> plotcommand.gnu
echo "set title '$EFFTITLE'" >> plotcommand.gnu
echo "set output '$EFFOUT'" >> plotcommand.gnu
echo "plot 'oldeffs.out' ls 3 title 'Initial Matrix-Matrix Multiply', 'neweffs.out' ls 2 title 'Matrix-Matrix Multiply With Optimzations From 3.1 and 3.2', 'new2effs.out' ls 4 title 'Matrix-Matrix Multiply With All Optimizations" >> plotcommand.gnu
fi

gnuplot plotcommand.gnu
rm -f plotcommand.gnu
if [ $1 == "time" ]
then
ps2pdf time.ps time.pdf
rm time.ps
else
ps2pdf eff.ps eff.pdf
rm eff.ps
fi
