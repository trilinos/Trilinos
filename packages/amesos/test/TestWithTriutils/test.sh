PROCS="1 2 4"
MATRICES="laplace_2d recirc_2d"
SIZE="100 10000 40000"

for p in $PROCS
do
  for m in $MATRICES
  do
    for s in $SIZE
    do
      mpirun -np $p TestOptions -problem_type=$m -problem_size=$s
    done
  done
done
