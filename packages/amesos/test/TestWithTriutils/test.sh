PROCS="1 2 4"
MATRICES="laplace_2d recirc_2d"
SIZE="100 10000 40000"

## Some machines use a command different than mpirun to run mpi jobs.  The
## test-harness.plx script sets the environment variable
## "TRILINOS_TEST_HARNESS_MPIGO_COMMAND".  We test for
## this value below.  If not set, we set it to a default value.

set mpigo = `printenv TRILINOS_TEST_HARNESS_MPIGO_COMMAND`

if ("$mpigo" == "") then
    set mpigo = "mpirun -np "
endif

for p in $PROCS
do
  for m in $MATRICES
  do
    for s in $SIZE
    do
      $mpigo $p TestOptions -problem_type=$m -problem_size=$s
    done
  done
done
