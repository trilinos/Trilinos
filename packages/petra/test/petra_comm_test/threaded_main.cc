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

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // I'm alive !!!



  Petra_Comm & petracomm = *new Petra_Comm( MPI_COMM_WORLD );
  int MyPID =  petracomm.MyPID();
  int NumProc =  petracomm.NumProc();
  cout << "Processor "<<MyPID<<" of " << NumProc << " is alive."<<endl;

  if (NumProc!=2) {
    cout << " This special routine only works for 2 processors " << endl;
    abort();
  }


  Petra_Comm * other_comm;
  MPI_Status status;

  unsigned int icomm = &petracomm;

  if (MyPID==1) cout << "Address of Petra_Comm object on PE 1 = " << &petracomm << endl;

  if (MyPID==0) MPI_Recv((void *) &other_comm, 1, MPI_UNSIGNED, 1, 99, MPI_COMM_WORLD, &status);
  else MPI_Send ( (void *) &icomm, 1, MPI_UNSIGNED, 0, 99, MPI_COMM_WORLD);

  if (MyPID==0) cout << "Address of other Petra_Comm object on PE 0 = " << other_comm << endl;
  
  int otherPID = other_comm->MyPID();

  if (MyPID==0) cout << "Processor "<<MyPID<<" of " << NumProc
		     << " has a neighbor processor with ID "
		     << otherPID << " of " << other_comm->NumProc() <<endl;
 
  delete &petracomm;
  MPI_Finalize();

  return 0;
}

