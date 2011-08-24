#ifndef SHYLU_H
#define SHYLU_H

#include "Epetra_CrsMatrix.h" 
#include "Epetra_Map.h" 
#include "Epetra_MultiVector.h" 
#include "Epetra_LinearProblem.h" 
#include "Epetra_SerialComm.h"
#include "Amesos_BaseSolver.h" 
#include "AztecOO.h"
#include "Isorropia_EpetraProber.hpp"

//#include "shylu_debug_manager.hpp"

#define MIN(a, b) (((a) < (b)) ? a : b)
#define MAX(a, b) (((a) > (b)) ? a : b)

// This is NOT just the symbolic structure, needs a better name
typedef struct
{
    Teuchos::RCP<Epetra_CrsMatrix> D;        // D Matrix
    Teuchos::RCP<Epetra_CrsMatrix> C;        // Column separator
    Teuchos::RCP<Epetra_CrsMatrix> R;        // Row separator
    Teuchos::RCP<Epetra_CrsMatrix> G;        // G Matrix (A22 block)
    Teuchos::RCP<Epetra_LinearProblem> LP;   // Local problem to solve D
    Teuchos::RCP<Amesos_BaseSolver> Solver;  // Local solver for D
    Teuchos::RCP<Epetra_CrsGraph> Sg;        // The approximate graph of S
                                             // Graph(S) + few diagonals
    Teuchos::RCP<Isorropia::Epetra::Prober> prober;  // Prober for Sbar
} shylu_symbolic;


typedef struct
{
    int Dnr;                    // #local rows
    int Dnc;                    // #local cols
    int Snr;                    // #remote rows
    int *DRowElems;             // local rows
    int *SRowElems;             // remote rows
    int *DColElems;             // Columns in D
    int *gvals;                 // O(n) array differentiating local/global
                                //  row/col
    //Epetra_SerialComm *SComm;   // Serial comm for block diagonals
    //Epetra_Map *LDRowMap;       // RowMap for block diagonals
    //Epetra_Map *LDColMap;       // ColMap for block diagonals
    //Epetra_CrsMatrix *D;        // Actual D Matrix, not reqd for Amesos_KLU
                                // but required for Amesos_Pardiso
    Teuchos::RCP<Epetra_CrsMatrix> Sbar; // Approx Schur complement
    AztecOO *innersolver;            // inner solver
    Epetra_LinearProblem *LP2;   // Local problem to solve
    Amesos_BaseSolver *dsolver;  // Local Subdomain solver
    int lmax;                    // May be this is optimizing too much
    int rmax;                    // May be this is optimizing too much
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
    int debug_level;
    //DebugManager dm;
} shylu_config;

int shylu_factor(Epetra_CrsMatrix *A, shylu_symbolic *ssym, shylu_data *data,
                shylu_config *config);

int shylu_symbolic_factor
(
    Epetra_CrsMatrix *A,    // i/p: A matrix
    shylu_symbolic *ssym,   // symbolic structure
    shylu_data *data,       // numeric structure, TODO: Required ?
    shylu_config *config   // i/p: library configuration
);

int shylu_solve(shylu_symbolic *ssym, shylu_data *data, shylu_config *config,
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
