// Petra_Object Test routine

#include <math.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#ifdef PETRA_MPI
#include <mpi.h>
#endif

#include "TSF_ParameterList.h"
//#include "Petra_Comm.h"

int main(int argc, char *argv[]) {
  /*
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
  */
  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


  // I'm alive !!!
  //if (verbose) cout << comm <<endl;

  TSF::ParameterList p;

  p.addParameter("Solver", "GMRES");
  p.addParameter("Preconditioner", "ILUT");
  p.addParameter("Threshold", 0.01); 
  p.addParameter("Max Fill", 20);
  p.addParameter("Scaling", "ColSum");

  char * solver;
  if (p.getParameter("Solver", solver))
    if (verbose) cout << "Solver selected is:  " << solver << endl;

  char * preconditioner;
  if (p.getParameter("Preconditioner", preconditioner)) {
    if (strcmp(preconditioner,"ILUT")==0) {
      if (verbose) cout << "Preconditioner selected is ILUT." << endl;
      double threshold = 0.1; // Default threshold
      if (p.getParameter("Threshold", threshold)) 
	if (verbose) 
	  cout << "Over-riding default threshold with threshold = " << threshold << endl;
      int maxfill = 10;
      if (p.getParameter("Max Fill", maxfill))
	if (verbose) 
	  cout << "Over-riding default maximum fill with maxfill = " << maxfill << endl;
    }
    else if (verbose) cout << "Preconditioner is not ILUT." << endl;
  }
  else if (verbose)
    cout << "No preconditioner selected.  Using default preconditioner." << endl;

  char * scaling;
  if (p.getParameter("Scaling", scaling)) {
    if (strcmp(scaling,"RowSum")==0) {
      if (verbose) cout << "Scaling selected is RowSum." << endl;
    }
    else if (verbose) cout << "Other scaling selected: " << scaling << endl;
  }
  else if (verbose) cout << "No scaling selected" << endl;
  
  
#ifdef PETRA_MPI
  MPI_Finalize();
#endif
  return 0;
}

/*
  end of file main.cc
*/
