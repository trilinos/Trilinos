#include "Teuchos_MPISession.hpp"
#include "Teuchos_TestForException.hpp"

using namespace Teuchos;

int MPISession::rank_ = 0 ;
int MPISession::nProc_ = 1 ;

void MPISession::init(int* argc, void*** argv)
{
#ifdef HAVE_MPI
	/* initialize MPI */
	int mpiHasBeenStarted = 0;
	MPI_Initialized(& mpiHasBeenStarted);
	int mpierr = 0 ;
	if (!mpiHasBeenStarted)
		{
			mpierr = ::MPI_Init (argc, (char ***) argv);
      TEST_FOR_EXCEPTION(mpierr != 0, runtime_error,
                         "Error code=" << mpierr 
                         << " detected in MPI_Init()");
		}
	
	/* find rank */
	mpierr = ::MPI_Comm_rank (MPI_COMM_WORLD, &rank_);
	TEST_FOR_EXCEPTION(mpierr != 0, runtime_error,
                     "Error code=" << mpierr 
                     << " detected in MPI_Comm_rank()");

	/* find number of procs */
	mpierr = ::MPI_Comm_size (MPI_COMM_WORLD, &nProc_);

	TEST_FOR_EXCEPTION(mpierr != 0, runtime_error,
                     "Error code=" << mpierr 
                     << " detected in MPI_Comm_size()");

  /* get machine name */
  int nameLen;
	char procName[MPI_MAX_PROCESSOR_NAME];
  mpierr = ::MPI_Get_processor_name(procName,&nameLen);

  TEST_FOR_EXCEPTION(mpierr != 0, runtime_error,
                     "Error code=" << mpierr 
                     << " detected in MPI_Get_processor_name()");

  cerr << "Teuchos::MPISession::init() started processor " << procName << endl;
  
#else
  cerr << "Teuchos::MPISession::init() started serial run" << endl;
#endif
}

void MPISession::finalize()
{
#ifdef HAVE_MPI
	int mpierr = ::MPI_Finalize();

	TEST_FOR_EXCEPTION(mpierr != 0, runtime_error,
                     "Error code=" << mpierr << " detected in MPI_Finalize()");
#endif
}
