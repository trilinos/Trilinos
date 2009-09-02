
/*TEST
PATH='tests/simplePi.c'
CCFLAGS=""
INPUT="10\n0\n"
OUTPUT='Enter the number of intervals: (0 quits)\npi is approximately 3.1424, Error is 0.0008\nEnter the number of intervals: (0 quits)\n'
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********       simplePi.c          ********************/
/******************************************************************/
/* Authors: William Gropp, Ewing Lusk                             */
/* Argonne National Laboratory                                    */
/* http://www-unix.mcs.anl.gov/mpi/tutorial/mpiintro/ppframe.htm  */
/******************************************************************/
/******************************************************************/

#include <math.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
  int done = 0, n, myid, numprocs, i, rc;
  double PI25DT = 3.141592653589793238462643;
  double mypi, pi, h, sum, x, a;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  while(!done) {
    if (myid == 0) {
      printf("Enter the number of intervals: (0 quits)\n");
      scanf("%d", &n);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (n == 0) break;
    h = 1.0/ (double) n;
    sum = 0.0;
    for (i = myid+1; i<=n; i+=numprocs) {
      x = h * ((double)i-0.5);
      sum += 4.0 / (1.0 + x*x);
    }
    mypi = h * sum;
    MPI_Reduce (&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid == 0)
      printf("pi is approximately %6.4f, Error is %6.4f\n", pi, fabs(pi-PI25DT));
      /* abs should be fabs */
  }
  MPI_Finalize();
  return 0;
}

    
