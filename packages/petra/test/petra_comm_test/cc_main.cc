// Petra_Comm Test routine

#include <math.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#ifdef PETRA_MPI
#include <mpi.h>
#endif
#include "Petra_Comm.h"
#include "Petra_Time.h"

int main(int argc, char *argv[]) {

#ifdef PETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


#ifdef PETRA_MPI

  Petra_Comm & petracomm = *new Petra_Comm( MPI_COMM_WORLD );
  int MyPID =  petracomm.MyPID();
  int NumProc =  petracomm.NumProc();

  assert(petracomm.MyPID()==rank);
  assert(petracomm.NumProc()==size);


  MPI_Comm MPIComm1 = petracomm.Comm();
  int size1, rank1;
  MPI_Comm_size(MPIComm1, &size1);
  MPI_Comm_rank(MPIComm1, &rank1);
  if (verbose) cout << petracomm <<  ".  Using MPI_Comm from Petra_Comm: "
                    << "Processor "<< rank1 <<" of " << size1
		          << " (should be the same)."<<endl;

  assert(rank1==rank);
  assert(size1==size);



  // Do some timing to test barrier function
  
  Petra_Time & before_barrier = *new Petra_Time(petracomm);
  Petra_Time & after_barrier = *new Petra_Time(petracomm);
  // Give each processor rank+1 amount of work
  // Time before barrier should increase roughly linearly
  // Time after barrier should be same for all processors
  double sum = 0.0;
  for (int j=0; j<rank+1; j++)
    for (int i=0; i<1000000; i++) sum += drand48();
  sum /= rank+1;
  if (verbose) cout << "Processor "<<MyPID
		    <<" Time to reach barrier: "
		    << before_barrier.ElapsedTime() << endl;
  petracomm.Barrier();
  if (verbose) cout << "Processor "<<MyPID << " Sum result  "
		    << sum <<" Time to get beyond barrier: "
		    << after_barrier.ElapsedTime() << endl;
 
  delete &before_barrier;
  delete &after_barrier;
  delete &petracomm;
  MPI_Finalize();

#else

  // Test serial interface first
  Petra_Comm & comm = *new Petra_Comm();
  int MyPID = comm.MyPID();
  int NumProc = comm.NumProc();
  if (verbose) cout << comm << endl;

  assert(MyPID==0);
  assert(NumProc==1);
  comm.Barrier();
  if (verbose) cout << comm << " is past serial barrier." << endl;
  delete &comm;

#endif

  return 0;
}

/*
  end of file main.cc
*/
