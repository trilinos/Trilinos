
#include "Epetra_CrsMatrix.h" 
#include "Epetra_MultiVector.h" 
#include "Epetra_LinearProblem.h" 
#include "Amesos_BaseSolver.h" 

int HyperLU_factor(Epetra_CrsMatrix *A, int sym, Epetra_MultiVector *&localS,
        Epetra_LinearProblem *&LP, Amesos_BaseSolver *&Solver, int &Dnr, 
        int *&DRowElems, int &Snr, int *&SRowElems, int *&piv,
        Epetra_MultiVector *&CMV);

