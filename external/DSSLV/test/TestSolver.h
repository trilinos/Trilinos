#include "Epetra_Comm.h"
enum SparseSolverType { UMFPACK, Aztec, SuperLU, SuperLUdist, SPOOLES, SPOOLESSERIAL, KUNDERT } ; 
enum DSS_MatrixType { DSS_Serial, DSS_Distributed } ; 

void TestSolver( Epetra_Comm &Comm, char *matrix_file, 
		 SparseSolverType SparseSolver,
		 bool tranpose, int special, DSS_MatrixType MatrixType );

void TestMultiSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special );

void TestMrhsSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special );


#if ( ! defined( CYGWINGCC ) && ! defined( TFLop ) ) 
#define TEST_KUNDERT
#endif
