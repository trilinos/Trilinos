
#include "Epetra_CrsMatrix.h" 
#include "Epetra_MultiVector.h" 
#include "Epetra_LinearProblem.h" 
#include "Amesos_BaseSolver.h" 

#define MIN(a, b) (((a) < (b)) ? a : b)
#define MAX(a, b) (((a) > (b)) ? a : b)

int HyperLU_factor(Epetra_CrsMatrix *A, int sym,
        Epetra_LinearProblem *&LP, Amesos_BaseSolver *&Solver, 
        Epetra_CrsMatrix *&Cptr, int &Dnr, 
        int *&DRowElems, int &Snr, int *&SRowElems,
        Teuchos::RCP<Epetra_CrsMatrix>& Sbar, double Sdiagfactor);

