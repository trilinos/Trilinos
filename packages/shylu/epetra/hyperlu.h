
#include "Epetra_CrsMatrix.h" 
#include "Epetra_MultiVector.h" 
#include "Epetra_LinearProblem.h" 
#include "Amesos_BaseSolver.h" 
#include "AztecOO.h"

#define MIN(a, b) (((a) < (b)) ? a : b)
#define MAX(a, b) (((a) > (b)) ? a : b)


typedef struct
{
    Epetra_LinearProblem *LP;   // Local problem to solve
    Amesos_BaseSolver *Solver;  // Local Subdomain solver
    Epetra_CrsMatrix *Cptr;     // Column separator
    int Dnr;                    // #local rows
    int Snr;                    // #remote rows
    int *DRowElems;             // local rows
    int *SRowElems;             // remote rows
    Teuchos::RCP<Epetra_CrsMatrix> Sbar; // Approx Schur complement
    AztecOO *innersolver;            // inner solver
} hyperlu_data;

typedef struct
{
    int sym;                    // flag for symmetry
    double Sdiagfactor;         // % of diagonals added to Schur complement
    int schurApproxMethod;      // ==1 implies blockdiagonal + A22
                                // ==2 implies dropping based

    double relative_threshold;  // Relative threshold for dropping
                                // only used if schurApproxMethod == 2

    int inner_maxiters;         // maximum iterations for inner solver
    double inner_tolerance;     // relative residual tolerance for inner solver
    string libName;             // library for the outer solver
} hyperlu_config;

int HyperLU_factor(Epetra_CrsMatrix *A, hyperlu_data *data, hyperlu_config 
                *config);

int hyperlu_solve(hyperlu_data *data, hyperlu_config *config,
    const Epetra_MultiVector& X, Epetra_MultiVector& Y);

Teuchos::RCP<Epetra_CrsMatrix> computeApproxSchur(hyperlu_config *config,
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap);
