// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef ML_NOX_PRECONDITIONER_UTILS_H
#define ML_NOX_PRECONDITIONER_UTILS_H

#include "Epetra_CrsGraph.h"
// the graph coloring class 
#include "Epetra_MapColoring.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#include "Epetra_CrsMatrix.h"


namespace ML_NOX
{
  Epetra_CrsGraph* deepcopy_graph(const Epetra_CrsGraph* oldgraph);
  Epetra_MapColoring* ML_Nox_collapsedcoloring(Epetra_CrsGraph* cgraph,
                                                     const int bsize,
                                                     bool diagonalonly,
						     int printlevel);
  Epetra_MapColoring* ML_Nox_standardcoloring(Epetra_CrsGraph* graph,
                                                      bool diagonalonly);
  bool Print_Epetra_CrsMatrix(Epetra_CrsMatrix& matrix);                                                    
  Epetra_CrsMatrix* StripZeros(Epetra_CrsMatrix& in, double eps = 0.0);                                                    
}




#endif
#endif
