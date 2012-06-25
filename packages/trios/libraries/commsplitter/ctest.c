#include <stdio.h>
#include <unistd.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    int np;
    int rank;
    char buf[256];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    gethostname(buf, 256);
    printf("ctest: rank %d of %d executing on %s\n", rank, np, buf);
    MPI_Finalize();

    return(0);
}
