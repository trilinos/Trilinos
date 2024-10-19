// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLU_SYMBOLIC_H
#define SHYLU_SYMBOLIC_H

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

#include "Epetra_CrsMatrix.h" 
#include "Epetra_Map.h" 
#include "Epetra_MultiVector.h" 
#include "Epetra_LinearProblem.h" 
#include "Isorropia_EpetraProber.hpp"
//#include "EpetraExt_Transpose_RowMatrix.h"
#include <EpetraExt_Reindex_LinearProblem2.h>
#include "Ifpack_Preconditioner.h"

// This is NOT just the symbolic structure, needs a better name
typedef struct
{
    Teuchos::RCP<Epetra_CrsMatrix> D;        // D Matrix
    //Teuchos::RCP<Epetra_CrsMatrix> DT;       // D Transpose Matrix
    //Teuchos::RCP<EpetraExt::RowMatrix_Transpose> transposer;
    Teuchos::RCP<Epetra_CrsMatrix> C;        // Column separator
    Teuchos::RCP<Epetra_CrsMatrix> R;        // Row separator
    Teuchos::RCP<Epetra_CrsMatrix> G;        // G Matrix (A22 block)
    Teuchos::RCP<Epetra_LinearProblem> LP;   // Local problem to solve D
    Teuchos::RCP<Epetra_LinearProblem> OrigLP;   // Local problem to solve D
    Teuchos::RCP<EpetraExt::ViewTransform<Epetra_LinearProblem> > ReIdx_LP ;
    Teuchos::RCP<Epetra_MultiVector> Dlhs;
    Teuchos::RCP<Epetra_MultiVector> Drhs;
    Teuchos::RCP<Epetra_MultiVector> Gvec;
    Teuchos::RCP<Amesos_BaseSolver> Solver;  // Local solver for D
    Teuchos::RCP<Ifpack_Preconditioner> ifSolver; //Local incomplete preconditioner
    Teuchos::RCP<Epetra_CrsGraph> Sg;        // The approximate graph of S
                                             // Graph(S) + few diagonals
    Teuchos::RCP<Isorropia::Epetra::Prober> prober;  // Prober for Sbar
} shylu_symbolic;


#endif // SHYLU_SYMBOLIC_H
