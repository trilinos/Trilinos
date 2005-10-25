
#include "ml_common.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef ML_NOX_PRECONDITIONER_UTILS_H
#define ML_NOX_PRECONDITIONER_UTILS_H

#include "Epetra_CrsGraph.h"
// the graph coloring class 
#include "Epetra_MapColoring.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"


namespace ML_NOX
{
  Epetra_CrsGraph* ML_NOX::deepcopy_graph(const Epetra_CrsGraph* oldgraph);
  Epetra_MapColoring* ML_NOX::ML_Nox_collapsedcoloring(Epetra_CrsGraph* cgraph,
                                                     Epetra_CrsGraph* graph, 
                                                     const int bsize,
                                                     bool diagonalonly,
						     int printlevel);
  Epetra_MapColoring* ML_NOX::ML_Nox_standardcoloring(Epetra_CrsGraph* graph,
                                                      bool diagonalonly);
}




#endif
#endif
