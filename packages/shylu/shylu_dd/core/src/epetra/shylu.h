// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file shylu.h

    \brief Main header file of ShyLU (Include main user calls)

    \author Siva Rajamanickam
*/
#ifndef SHYLU_H
#define SHYLU_H

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

// Epetra include
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_SerialComm.h"

// Amesos includes
#include "Amesos_BaseSolver.h"

// Ifpack includes
#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"

// AztecOO includes
#include "AztecOO.h"

// Isorropia includes
#include "Isorropia_EpetraProber.hpp"

// Amesos2 includes
#ifdef HAVE_SHYLU_DDCORE_AMESOS2
#include <Amesos2.hpp>
#endif

// Tpetra includes
#ifdef HAVE_SHYLU_DDCORE_TPETRA
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#endif

// Zoltan2 includes
#if defined(HAVE_SHYLU_DDCORE_ZOLTAN2CORE)
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#endif


// Shylu includes
#include "shylu_symbolic.h"
#include "shylu_config.h"
#include "shylu_probing_operator.h"
#include "shylu_amesos_schur_operator.h"

#include <IQRSolver.h>


#define SHYLU_CORE_MIN(a, b) (((a) < (b)) ? a : b)
#define SHYLU_CORE_MAX(a, b) (((a) > (b)) ? a : b)

/** \brief Main data structure holding needed offset and temp variables
 *
 * This structur contains ...
 */
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

/** \brief Main function call into ShylU
 *
 * How to use?
 */
int shylu_factor(Epetra_CrsMatrix *A, shylu_symbolic *ssym, shylu_data *data,
                shylu_config *config);

/** \brief Call symbolic factorization on matrix
 *
 */
int shylu_symbolic_factor
(
    Epetra_CrsMatrix *A,    // i/p: A matrix
    shylu_symbolic *ssym,   // symbolic structure
    shylu_data *data,       // numeric structure, TODO: Required ?
    shylu_config *config   // i/p: library configuration
);

/** \brief Call solve on multiple RHS
 *
 */
int shylu_solve(shylu_symbolic *ssym, shylu_data *data, shylu_config *config,
    const Epetra_MultiVector& X, Epetra_MultiVector& Y);

/** \brief Compute an approximate Schur Complement (Narrow Sep)
 *
 *  Computate an approximate Schur Complement either using ...
 */

Teuchos::RCP<Epetra_CrsMatrix> computeApproxSchur(shylu_config *config,
    shylu_symbolic *ssym,
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver,
    Ifpack_Preconditioner *ifSolver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap);

/** \brief Compute an approximate Shur Complete (Wide Sep)
 *
 * Compute an approximate Schur Complement based on a wide seperator.
 * Options include ...
 */
Teuchos::RCP<Epetra_CrsMatrix> computeApproxWideSchur(
    shylu_config *config,
    shylu_symbolic *ssym,   // symbolic structure
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver,
    Ifpack_Preconditioner *ifSolver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap);

/** \brief Compute an approximate Schur Complement using the option of Guided Probing
 *
 *  Compute an approximate Schur Complement based on probing of important nonzero values.
 */
Teuchos::RCP<Epetra_CrsMatrix> computeSchur_GuidedProbing
(
    shylu_config *config,
    shylu_symbolic *ssym,   // symbolic structure
    shylu_data *data,       // numeric structure
    Epetra_Map *localDRowMap
);
#endif // SHYLU_H
