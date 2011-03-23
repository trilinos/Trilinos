/** \file hyperlu_util.h

    \brief Utilities for HyperLU

    \author Siva Rajamanickam

*/
#ifndef HYPERLU_UTIL_H
#define HYPERLU_UTIL_H

#include <assert.h>
#include <mpi.h>
#include <iostream>
#include <sstream>

// To dump all the matrices into files.
//#define DUMP_MATRICES

#include "Isorropia_config.h" // Just for HAVE_MPI

// Epetra includes
#ifdef HAVE_MPI
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

using namespace std;

Epetra_CrsMatrix* balanceAndRedistribute(Epetra_CrsMatrix* A, 
                        Teuchos::ParameterList isoList);

void checkMaps(Epetra_CrsMatrix *A);

void findLocalColumns(Epetra_CrsMatrix *A, int *gvals, int &SNumGlobalCols);

void findNarrowSeparator(Epetra_CrsMatrix *A, int *gvals);

void findBlockElems(int nrows, int *rows, int *gvals, int Lnr, int *LeftElems, 
        int Rnr, int *RightElems, string s1, string s2);

#endif //HYPERLU_UTIL_H
