// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file shylu_util.h

    \brief Utilities for ShyLU

    \author Siva Rajamanickam

*/
#ifndef SHYLU_UTIL_H
#define SHYLU_UTIL_H

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

#include <assert.h>
#include <mpi.h>
#include <iostream>
#include <sstream>

// To dump all the matrices into files.
//#define DUMP_MATRICES

//#define TIMING_OUTPUT
//
#include "ShyLU_DDCore_config.h"

// Epetra includes
#ifdef HAVE_SHYLU_DDCORE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_Time.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_CrsMatrixIn.h"

// Amesos includes
#include "Amesos.h"
#include "Amesos_BaseSolver.h"

// AztecOO includes
#include "AztecOO.h"

#include "Isorropia_EpetraProber.hpp"


Epetra_CrsMatrix* balanceAndRedistribute(Epetra_CrsMatrix* A,
                        Teuchos::ParameterList isoList);

void checkMaps(Epetra_CrsMatrix *A);

void findLocalColumns(Epetra_CrsMatrix *A, int *gvals, int &SNumGlobalCols);

void findNarrowSeparator(Epetra_CrsMatrix *A, int *gvals);

void findBlockElems(Epetra_CrsMatrix *A, int nrows, int *rows, int *gvals,
        int Lnr, int *LeftElems,
        int Rnr, int *RightElems, std::string s1, std::string s2, bool cols);

#ifdef SHYLU_DEBUG

#define SHYLU_CORE_ASSERT(A) assert(A)

#else

#define SHYLU_CORE_ASSERT(A)

#endif

#endif //SHYLU_UTIL_H
