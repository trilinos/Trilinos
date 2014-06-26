#!/bin/bash


rm -f example*.txt
rm -f aggs*.txt # clean up
rm -f nodes*.txt

function genTitle() {
  title=" === $EXAMPLE ($NX x $NY) === \n\n solver: $XMLFILE   procs: $PROCS   Multigrid sweeps: $SWEEPS\n\n"
  echo -e $title
}

function printSubMenu() {
  clear
  title="Choose action"
  prompt="Select:"
  options=("rerun example" "change multigrid sweeps" "change mesh" "change solver" "change procs" "show output" "plot exact solution" "plot Multigrid solution" "plot error" "plot proc distribution")

  EXAMPLE=$4
  PROCS=2
  SWEEPS=1
  NX=$1
  NY=$2
  XMLFILE= #tutorial2a.xml
  EXE=$3

  genTitle $EXAMPLE $NX $NY $XMLFILE $PROCS $SWEEPS

  #echo "$title"
  PS3="$prompt "
  select opt in "${options[@]}" "Quit"; do

      case "$REPLY" in

      1 ) # clean up
          echo "No of multigrid sweeps=$SWEEPS"
	  rm -f example*.txt # clean up
	  rm -f aggs*.txt # clean up
	  rm -f nodes*.txt
	  rm -f output.log
	  runExample $NX $NY $SWEEPS $XMLFILE $PROCS $EXE
	  cat example*.txt > example.txt
	  clear
	  genTitle $EXAMPLE $NX $NY $XMLFILE $PROCS $SWEEPS
	  i=1
	  for value in "${options[@]}" "Quit"; do
	    printf "%-25s\n" "$i) $value"
	    i=$((i+1))
	  done | column
          ;;
      2 ) #clean up
	  rm -f example*.txt # clean up
	  rm -f aggs*.txt # clean up
	  rm -f nodes*.txt
	  rm -f output.log
	  echo
	  echo -n "Number of multigrid sweeps: "
	  read SWEEPS
	  echo
	  runExample $NX $NY $SWEEPS $XMLFILE $PROCS $EXE
	  cat example*.txt > example.txt
	  clear
	  genTitle $EXAMPLE $NX $NY $XMLFILE $PROCS $SWEEPS
	  i=1
	  for value in "${options[@]}" "Quit"; do
	    printf "%-25s\n" "$i) $value"
	    i=$((i+1))
	  done | column
          ;;

      3 ) #change mesh
      	  rm -f example*.txt # clean up
      	  rm -f aggs*.txt # clean up
      	  rm -f nodes*.txt
	  rm -f output.log
	  echo
          echo "*** Change mesh ***"
          echo ""
          echo -n "nx="
          read NX
          echo -n "ny="
          read NY
          echo ""
          echo -n "Prepare and solve example"
	  runExample $NX $NY $SWEEPS $XMLFILE $PROCS $EXE
	  echo " DONE"
	  echo
	  cat example*.txt > example.txt
  	  clear
	  genTitle $EXAMPLE $NX $NY $XMLFILE $PROCS $SWEEPS
	  i=1
	  for value in "${options[@]}" "Quit"; do
	    printf "%-25s\n" "$i) $value"
	    i=$((i+1))
	  done | column
          ;;
      4 ) #change solver
	  rm -f example*.txt # clean up
	  rm -f aggs*.txt # clean up
	  rm -f nodes*.txt
	  rm -f output.log
	  echo
          echo "*** Change solver ***"
          echo ""
          echo -n "XML file="
          read XMLFILE
          echo ""
          echo -n "Prepare and solve example"
	  runExample $NX $NY $SWEEPS $XMLFILE $PROCS $EXE
	  echo " DONE"
	  echo
	  cat example*.txt > example.txt
  	  clear
	  genTitle $EXAMPLE $NX $NY $XMLFILE $PROCS $SWEEPS
	  i=1
	  for value in "${options[@]}" "Quit"; do
	    printf "%-25s\n" "$i) $value"
	    i=$((i+1))
	  done | column
          ;;
      5 ) #change procs
	  rm -f example*.txt # clean up
	  rm -f aggs*.txt # clean up
	  rm -f nodes*.txt
	  rm -f output.log
	  echo
          echo "*** Change processors ***"
          echo ""
          echo -n "#processors="
          read PROCS
          echo ""
          echo -n "Prepare and solve example"
	  runExample $NX $NY $SWEEPS $XMLFILE $PROCS $EXE
	  echo " DONE"
	  echo
	  cat example*.txt > example.txt
  	  clear
	  genTitle $EXAMPLE $NX $NY $XMLFILE $PROCS $SWEEPS
	  i=1
	  for value in "${options[@]}" "Quit"; do
	    printf "%-25s\n" "$i) $value"
	    i=$((i+1))
	  done | column
          ;;
      6 ) less output.log
          ;;
      7 ) plotSolution $NX $NY "example.txt" $SWEEPS $XMLFILE $PROCS &
	  ;;
      8 ) plotMultigridSolution $NX $NY "example.txt" $SWEEPS $XMLFILE $PROCS &
	  ;;
      9 ) plotError $NX $NY "example.txt" $SWEEPS $XMLFILE $PROCS &
	  ;;
      10 ) plotProc $NX $NY "example.txt" $SWEEPS $XMLFILE $PROCS &
          ;;
      $(( ${#options[@]}+1 )) )
	  echo "Clear example.";
	  rm -f example*.txt # clean up
	  rm -f aggs*.txt # clean up
	  rm -f nodes*.txt
	  rm -f output.log
	  exit
	  break;;
      *) echo "Invalid option. Try another one.";continue;;

      esac

  done
}

function printMenu() {
  title="Select example"
  prompt="Pick an option:"
  options=("2D Laplace (50x50 mesh)" "2D Laplace (interactive)" "2D Recirc (50x50 mesh)" "2D Recirc")

  echo "$title"
  PS3="$prompt "
  select opt in "${options[@]}" "Quit"; do

      case "$REPLY" in

      1 ) NX=50
          NY=50
	  runExample $NX $NY 1 xml/s2a.xml 2 ./MueLu_laplace2d.exe "2D Laplace"
	  cat example*.txt > example.txt
	  printSubMenu $NX $NY "./MueLu_laplace2d.exe" "2D Laplace"
          ;;
      2 ) echo ""
          echo "*** 2D Laplace example ***"
          echo ""
          echo -n "nx="
          read NX
          echo -n "ny="
          read NY
          echo ""
          echo -n "Prepare and solve example"
	  runExample $NX $NY 1 xml/s2a.xml 2 "./MueLu_laplace2d.exe" "2D Laplace"
	  echo " DONE"
	  echo
	  cat example*.txt > example.txt
	  printSubMenu $NX $NY "./MueLu_laplace2d.exe" "2D Laplace"
	  ;;
      3 ) echo ""
          echo "*** 2D Recirc example ***"
          echo ""
          NX=50
          NY=50
          echo ""
          echo -n "Prepare and solve example"
	  runExample $NX $NY 1 xml/s3a.xml 2 "./MueLu_recirc2d.exe" "2D Recirc"
	  echo " DONE"
	  echo
	  cat example*.txt > example.txt
	  printSubMenu $NX $NY "./MueLu_recirc2d.exe" "2D Recirc"
	  ;;
      4 ) echo ""
          echo "*** 2D Recirc example ***"
          echo ""
          echo -n "nx="
          read NX
          echo -n "ny="
          read NY
          echo ""
          echo -n "Prepare and solve example"
	  runExample $NX $NY 1 xml/s3a.xml 2 "./MueLu_recirc2d.exe" "2D Recirc"
	  echo " DONE"
	  echo
	  cat example*.txt > example.txt
	  printSubMenu $NX $NY "./MueLu_recirc2d.exe" "2D Recirc"
	  ;;
      $(( ${#options[@]}+1 )) )
	  echo "Goodbye!";
	  rm -f example*.txt # clean up
	  rm -f aggs*.txt # clean up
	  rm -f output.log
	  break;;
      *) echo "Invalid option. Try another one.";continue;;

      esac

  done
}

function runExample() {
  mpirun -np $5 $6 --nx=$1 --ny=$2 --mgridSweeps=$3 --xml=$4 | tee output.log
  echo "Press q to return." >> output.log
  python ./MueLu_Agg2VTK.py &
  rm -f nodes*.txt
}

function plotSolution() {
  `gnuplot -persist << _TTT_
  set dgrid3d $2,$1
  set style data lines
  set nolabel
  set key off
  set autoscale
  splot "$3" using 3:4:5
  quit
_TTT_`
#set title "Exact solution"
}

function plotMultigridSolution() {
  `gnuplot -persist << _TTT_
  set dgrid3d $2,$1
  set style data lines
  set nolabel
  set key off
  set autoscale
  splot "$3" using 3:4:7
  quit
_TTT_`
  #set title "Multigrid solution after $SWEEPS V-cycle sweeps"
}

function plotError() {
  #echo $VAR | gnuplot -persist
  `gnuplot -persist << _TTT_
  set dgrid3d $2,$1
  set style data lines
  set palette model RGB defined ( 0 'red', 1 'green', 2 'blue', 3 'yellow', 4 'pink')
  set nolabel
  set key off
  set autoscale
  set hidden3d
  splot "$3" using 3:4:($ 5-$ 7)
_TTT_`
  #set title "Error after $SWEEPS V-cycle sweeps"
}

function plotProc() {
  #echo $VAR | gnuplot -persist
  `gnuplot -persist << _TTT_
  set dgrid3d $2,$1
  set style data lines
  set palette model RGB defined ( 0 'red', 1 'green', 2 'blue', 3 'yellow', 4 'pink')
  set nolabel
  set key off
  set autoscale
  splot "$3" using 3:4:1 with points palette
_TTT_`
  #set title "Distribution of processors"
}

printMenu




