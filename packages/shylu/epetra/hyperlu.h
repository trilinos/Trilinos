
#include "Epetra_CrsMatrix.h" 
#include "Epetra_MultiVector.h" 
#include "Epetra_LinearProblem.h" 
#include "Amesos_BaseSolver.h" 

#define MIN(a, b) (((a) < (b)) ? a : b)
#define MAX(a, b) (((a) > (b)) ? a : b)


typedef struct
{
    Epetra_LinearProblem *LP;  // Local problem to solve
    Amesos_BaseSolver *Solver; // Local Subdomain solver
    Epetra_CrsMatrix *Cptr;    // Column separator
    int Dnr;                   // #local rows
    int Snr;                   // #remote rows
    int *DRowElems;             // local rows
    int *SRowElems;             // remote rows
    Teuchos::RCP<Epetra_CrsMatrix> Sbar; // Approx Schur complement
} hyperlu_data;

int HyperLU_factor(Epetra_CrsMatrix *A, int sym, hyperlu_data *data,
        double Sdiagfactor);
/*int HyperLU_factor(Epetra_CrsMatrix *A, int sym,
        Epetra_LinearProblem *&LP, Amesos_BaseSolver *&Solver, 
        Epetra_CrsMatrix *&Cptr, int &Dnr, 
        int *&DRowElems, int &Snr, int *&SRowElems,
        Teuchos::RCP<Epetra_CrsMatrix>& Sbar, double Sdiagfactor);*/
