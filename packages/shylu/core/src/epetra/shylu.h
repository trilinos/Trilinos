
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef SHYLU_H
#define SHYLU_H

#include "Epetra_CrsMatrix.h" 
#include "Epetra_Map.h" 
#include "Epetra_MultiVector.h" 
#include "Epetra_LinearProblem.h" 
#include "Epetra_SerialComm.h"
#include "Amesos_BaseSolver.h"
#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "AztecOO.h"
#include "Isorropia_EpetraProber.hpp"

#include "shylu_symbolic.h"
#include "shylu_config.h"
#include "shylu_probing_operator.h"
#include "shylu_amesos_schur_operator.h"

#include <IQRSolver.h>

//#include "shylu_debug_manager.hpp"

#define MIN(a, b) (((a) < (b)) ? a : b)
#define MAX(a, b) (((a) > (b)) ? a : b)

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
    Teuchos::RCP<Epetra_Map> LDRowMap;       // RowMap for block diagonals
    Teuchos::RCP<Epetra_Map> LGRowMap;       // RowMap for G (local)
    Teuchos::RCP<Epetra_Map> GMap;           // Dist Map for G

    Teuchos::RCP<Epetra_Import> BdImporter;
    Teuchos::RCP<Epetra_Import> DistImporter;
    Teuchos::RCP<Epetra_Import> BsImporter;
    Teuchos::RCP<Epetra_Import> XsImporter;
    Teuchos::RCP<Epetra_Export> XdExporter;
    Teuchos::RCP<Epetra_Export> XsExporter;

    Teuchos::RCP<Epetra_MultiVector> localrhs;
    Teuchos::RCP<Epetra_MultiVector> temp1;
    Teuchos::RCP<Epetra_MultiVector> temp2;
    Teuchos::RCP<Epetra_MultiVector> Bs;
    Teuchos::RCP<Epetra_MultiVector> Xs;
    Teuchos::RCP<Epetra_MultiVector> LocalXs;
    Teuchos::RCP<Epetra_MultiVector> temp3;
    Teuchos::RCP<Epetra_MultiVector> locallhs;

    // temp timers
    //Teuchos::RCP<Teuchos::Time> importExportTime;
    //Teuchos::RCP<Teuchos::Time> innerIterTime;
    //Teuchos::RCP<Teuchos::Time> fwdTime;
    //Teuchos::RCP<Teuchos::Time> amesosSchurTime;

    //Epetra_CrsMatrix *D;        // Actual D Matrix, not reqd for Amesos_KLU
                                // but required for Amesos_Pardiso
    Teuchos::RCP<IQR::IQRSolver> iqrSolver; // Solver object for IQR func.
    Teuchos::RCP<Epetra_CrsMatrix> Sbar; // Approx Schur complement
    Teuchos::RCP<Epetra_CrsGraph> localSbargraph; // graph of local Sbar
    AztecOO *innersolver;            // inner solver
    Teuchos::RCP<Epetra_MultiVector> Sbarlhs;
    Teuchos::RCP<Epetra_MultiVector> Sbarrhs;
    Teuchos::RCP<Epetra_LinearProblem> LP2;   // Local problem to solve
    Teuchos::RCP<Epetra_LinearProblem> OrigLP2;   // Local problem to solve D
	Teuchos::RCP<EpetraExt::ViewTransform<Epetra_LinearProblem> > ReIdx_LP2;
	Amesos_BaseSolver *dsolver;  // Local Subdomain solver
    Teuchos::RCP<Ifpack_Preconditioner> schur_prec;
    Teuchos::RCP<ShyLU_Probing_Operator> schur_op;
    int lmax;                    // May be this is optimizing too much
    int rmax;                    // May be this is optimizing too much
    Teuchos::RCP<Isorropia::Epetra::Prober> guided_prober;  // Guided Prober for Sbar
    int num_compute;            // # of times Compute() has been called before
                                // or in otherwords #nonlinear iteration-1
} shylu_data;

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
    shylu_symbolic *ssym,
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver,
    Ifpack_Preconditioner *ifSolver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap);

Teuchos::RCP<Epetra_CrsMatrix> computeApproxWideSchur(
    shylu_config *config,
    shylu_symbolic *ssym,   // symbolic structure
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver,
    Ifpack_Preconditioner *ifSolver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap);

Teuchos::RCP<Epetra_CrsMatrix> computeSchur_GuidedProbing
(
    shylu_config *config,
    shylu_symbolic *ssym,   // symbolic structure
    shylu_data *data,       // numeric structure
    Epetra_Map *localDRowMap
);
#endif // SHYLU_H
