#include "Epetra_Comm.h"
enum SparseSolverType { UMFPACK, Aztec, SuperLU, SuperLUdist, SPOOLES, SPOOLESSERIAL, KUNDERT } ; 
enum AMESOS_MatrixType { AMESOS_Serial, AMESOS_Distributed } ; 

int TestSolver( Epetra_Comm &Comm, char *matrix_file, 
		 SparseSolverType SparseSolver,
		 bool tranpose, int special, AMESOS_MatrixType MatrixType );

int TestMultiSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special );

int TestMrhsSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special );

#if 0
#if ( ! defined( CYGWINGCC ) && ! defined( TFLOP ) ) 
#define TEST_KUNDERT
#endif
#if ( ! defined( TFLOP ) ) 
#define TEST_AZTEC
#endif
#endif
