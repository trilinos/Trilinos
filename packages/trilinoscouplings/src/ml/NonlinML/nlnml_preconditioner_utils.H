// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file nlnml_preconditioner_utils.H
 *
 * \brief ML nonlinear preconditioner and solver utility routines
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef NLNML_PRECONDITIONER_UTILS_H
#define NLNML_PRECONDITIONER_UTILS_H

#include "Epetra_CrsGraph.h"
// the graph coloring class 
#include "Epetra_MapColoring.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RefCountPtr.hpp"

using namespace Teuchos;
namespace NLNML
{
  Epetra_CrsGraph* deepcopy_graph(const Epetra_CrsGraph* oldgraph);
  RefCountPtr<Epetra_MapColoring> Collapsedcoloring(Epetra_CrsGraph* cgraph,
                                                    const int bsize,
                                                    bool diagonalonly,
			 		            int printlevel);
  RefCountPtr<Epetra_MapColoring> Standardcoloring(Epetra_CrsGraph* graph,
                                                   bool diagonalonly);
  bool Print_Epetra_CrsMatrix(Epetra_CrsMatrix& matrix);                                                    
  RefCountPtr<Epetra_CrsMatrix> StripZeros(RefCountPtr<Epetra_CrsMatrix> in, double eps = 0.0);                                                    
}




#endif
#endif
