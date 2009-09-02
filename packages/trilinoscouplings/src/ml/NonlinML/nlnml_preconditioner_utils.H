/*
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
*/
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
