#include "Epetra_Comm.h"
enum SparseSolverType { UMFPACK, Aztec, SuperLU, SuperLUdist, SPOOLES, SPOOLESSERIAL, KUNDERT } ; 
enum AMESOS_MatrixType { AMESOS_Serial, AMESOS_Distributed } ; 

int TestSolver( Epetra_Comm &Comm, char *matrix_file, 
		 SparseSolverType SparseSolver,
		 bool tranpose, int special, AMESOS_MatrixType MatrixType );

void TestMe( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special );

void TestMultiSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special );

void TestMrhsSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special );


#if ( ! defined( CYGWINGCC ) && ! defined( TFLop ) ) 
#define TEST_KUNDERT
#endif
