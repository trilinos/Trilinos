// Petra_Comm Test routine

#include <math.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <iostream.h>
#include <strstream.h>
#include <stdio.h>
#include <mpi.h>
#include "az_aztec.h"

extern "C" void AZ_sync1(int proc_config[]);

int main(int argc, char *argv[]) {


  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int    proc_config[AZ_PROC_SIZE];/* Processor information. */
  AZ_set_proc_config(proc_config,MPI_COMM_WORLD);
  int MyPID = proc_config[AZ_node];
  int NumProc = proc_config[AZ_N_procs];


  bool verbose = true;


  // Do some timing to test barrier function
  // Give each processor rank+1 amount of work
  // Time before barrier should increase roughly linearly
  // Time after barrier should be about the same for all processors
  double starttime = MPI_Wtime();
  double sum = 0.0;
  for (int j=0; j<rank+1; j++)
    for (int i=0; i<1000000; i++) sum += drand48();
  sum /= rank+1;
  if (verbose) cout << "Processor "<<MyPID
		    <<" Time to reach barrier: "
		    << MPI_Wtime() - starttime << endl;
  AZ_sync1(proc_config);
  if (verbose) cout << "Processor "<<MyPID << " Sum result  "
		    << sum <<" Time to get beyond barrier: "
		    << MPI_Wtime() - starttime << endl;
 
  MPI_Finalize();

  return 0;
}

