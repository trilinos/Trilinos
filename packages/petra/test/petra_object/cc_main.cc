// Petra_Object Test routine

#include <math.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#ifdef PETRA_MPI
#include <mpi.h>
#endif

#include "Petra_Object.h"
#include "Petra_Comm.h"

int main(int argc, char *argv[]) {

#ifdef PETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Petra_Comm comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

  Petra_Comm comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


  // I'm alive !!!
  if (verbose) cout << comm <<endl;

  Petra_Object obj;

  if (verbose) cout << "This is the default Petra_Object Name: " << obj <<endl;

  obj.SetLabel("New name for Petra_Object");

  if (verbose) cout << "This should say \"New name for Petra_Object\": " << obj <<endl;

#ifdef PETRA_MPI
  MPI_Finalize();
#endif
  return 0;
}

/*
  end of file main.cc
*/
