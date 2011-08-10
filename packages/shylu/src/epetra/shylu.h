#ifndef SHYLU_H
#define SHYLU_H

#include "Epetra_CrsMatrix.h" 
#include "Epetra_Map.h" 
#include "Epetra_MultiVector.h" 
#include "Epetra_LinearProblem.h" 
#include "Epetra_SerialComm.h"
#include "Amesos_BaseSolver.h" 
#include "AztecOO.h"

#define MIN(a, b) (((a) < (b)) ? a : b)
#define MAX(a, b) (((a) > (b)) ? a : b)


typedef struct
{
    Epetra_LinearProblem *LP;   // Local problem to solve
    Amesos_BaseSolver *Solver;  // Local Subdomain solver
    Epetra_CrsMatrix *Cptr;     // Column separator
    Epetra_CrsMatrix *Rptr;     // Row separator
    int Dnr;                    // #local rows
    int Snr;                    // #remote rows
    int *DRowElems;             // local rows
    int *SRowElems;             // remote rows
    //Epetra_SerialComm *SComm;   // Serial comm for block diagonals
    //Epetra_Map *LDRowMap;       // RowMap for block diagonals
    //Epetra_Map *LDColMap;       // ColMap for block diagonals
    //Epetra_CrsMatrix *D;        // Actual D Matrix, not reqd for Amesos_KLU
                                // but required for Amesos_Pardiso
    Teuchos::RCP<Epetra_CrsMatrix> Sbar; // Approx Schur complement
    AztecOO *innersolver;            // inner solver
    Epetra_LinearProblem *LP2;   // Local problem to solve
    Amesos_BaseSolver *dsolver;  // Local Subdomain solver
} shylu_data;

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
    string schurSolver;         // Solver for the Schur complement
    int sep_type;
} shylu_config;

int shylu_factor(Epetra_CrsMatrix *A, shylu_data *data, shylu_config 
                *config);

int shylu_solve(shylu_data *data, shylu_config *config,
    const Epetra_MultiVector& X, Epetra_MultiVector& Y);

Teuchos::RCP<Epetra_CrsMatrix> computeApproxSchur(shylu_config *config,
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap);

Teuchos::RCP<Epetra_CrsMatrix> computeApproxWideSchur(shylu_config *config,
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap);

#endif // SHYLU_H
