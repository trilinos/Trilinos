#ifndef TEUCHOS_MPISESSION_H
#define TEUCHOS_MPISESSION_H

#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

namespace Teuchos
{
  /**
   * Class MPISession provides methods for initializing, finalizing, 
   * and querying the global MPI session. 
   */
  class MPISession
    {
    public:
      /** initializer, calls MPI_Init() if necessary */
      static void init(int* argc, void*** argv);

      /** returns the process rank relative to MPI_COMM_WORLD */
      static int getRank() {return rank_;}

      /** returns the number of processors in MPI_COMM_WORLD */
      static int getNProc() {return nProc_;}

      /** finalizer, calls MPI_Finalize() if necessary */
      static void finalize();
    private:
      static int rank_;
      static int nProc_;
    };
}

#endif
