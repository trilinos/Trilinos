#include "Epetra_Comm.h"
enum SparseSolverType { UMFPACKOLD, Aztec, SuperLU, SuperLUdist, 
			SuperLUdist2, DSCPACK, DSCPACKOLD, UMFPACK, 
			SPOOLES, SPOOLESSERIAL, KUNDERT, MUMPS, KLU,
                        SUPERLUDIST } ; 
enum AMESOS_MatrixType { AMESOS_Serial, AMESOS_Distributed } ; 

int Amesos_TestSolver( Epetra_Comm &Comm, char *matrix_file, 
		 SparseSolverType SparseSolver,
		 bool tranpose, int special, AMESOS_MatrixType MatrixType );

int Amesos_TestMultiSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special, AMESOS_MatrixType MatrixType );

int Amesos_TestMrhsSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special, AMESOS_MatrixType MatrixType );

