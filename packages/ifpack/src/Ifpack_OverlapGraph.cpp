/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_OverlapGraph.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"

#ifdef HAVE_IFPACK_TEUCHOS
#include <Teuchos_ParameterList.hpp>
#include <ifp_parameters.h>
#endif

//==============================================================================
Ifpack_OverlapGraph::Ifpack_OverlapGraph(const Epetra_CrsGraph * UserMatrixGraph, int OverlapLevel)
  : OverlapGraph_(0),
    UserMatrixGraph_(UserMatrixGraph),
    UserMatrix_(0),
    OverlapRowMap_(0),
    OverlapLevel_(OverlapLevel)
{
  // Test for non-trivial overlap here so we can use it later.
  IsOverlapped_ = (OverlapLevel>0 && UserMatrixGraph->DomainMap().DistributedGlobal());

  ConstructOverlapGraph(UserMatrixGraph);

}
//==============================================================================
Ifpack_OverlapGraph::Ifpack_OverlapGraph(const Epetra_RowMatrix * UserMatrix, int OverlapLevel)
  : OverlapGraph_(0),
    UserMatrixGraph_(0),
    UserMatrix_(UserMatrix),
    OverlapRowMap_(0),
    OverlapLevel_(OverlapLevel)
{
  // Test for non-trivial overlap here so we can use it later.
  IsOverlapped_ = (OverlapLevel>0 && UserMatrix->OperatorDomainMap().DistributedGlobal());
  
  throw ReportError("This constructor is not implemented yet.  Need to add Epetra_SrcObject support to Epetra_Import/Export", -1);
}
//==============================================================================
Ifpack_OverlapGraph::Ifpack_OverlapGraph(const Ifpack_OverlapGraph & Source)
  : OverlapGraph_(Source.OverlapGraph_),
    UserMatrixGraph_(Source.UserMatrixGraph_),
    UserMatrix_(Source.UserMatrix_),
    OverlapRowMap_(Source.OverlapRowMap_),
    OverlapLevel_(Source.OverlapLevel_),
    IsOverlapped_(Source.IsOverlapped_)
{
  if (IsOverlapped_) {
    if (OverlapGraph_!=0) OverlapGraph_ = new Epetra_CrsGraph(*OverlapGraph_);
    if (OverlapRowMap_!=0) OverlapRowMap_ = new Epetra_BlockMap(*OverlapRowMap_);
  }
  
}
//==============================================================================
Ifpack_OverlapGraph::~Ifpack_OverlapGraph() {
  
  if (IsOverlapped_) {
    if (OverlapGraph_!=0) delete OverlapGraph_;
    if (OverlapRowMap_!=0) delete OverlapRowMap_;
  }
}

#ifdef HAVE_IFPACK_TEUCHOS
//==========================================================================
int Ifpack_OverlapGraph::SetParameters(const Teuchos::ParameterList& parameterlist,
				       bool cerr_warning_if_unused)
{
  Ifpack::param_struct params;
  params.int_params[Ifpack::level_overlap-FIRST_INT_PARAM] = OverlapLevel_;

  Ifpack::set_parameters(parameterlist, params, cerr_warning_if_unused);

  OverlapLevel_ = params.int_params[Ifpack::level_overlap-FIRST_INT_PARAM];
  return(0);
}
#endif

//==============================================================================
int Ifpack_OverlapGraph::ConstructOverlapGraph(const Epetra_CrsGraph * UserMatrixGraph) {

  OverlapGraph_ = (Epetra_CrsGraph *) UserMatrixGraph;
  OverlapRowMap_ = (Epetra_BlockMap *) &UserMatrixGraph->RowMap();

  if (!IsOverlapped_) return(0); // All done

  Epetra_CrsGraph * OldGraph;
  Epetra_BlockMap * OldRowMap;
  Epetra_BlockMap * DomainMap = (Epetra_BlockMap *) &UserMatrixGraph->DomainMap();
  Epetra_BlockMap * RangeMap = (Epetra_BlockMap *) &UserMatrixGraph->RangeMap();

  for (int level=1; level <= OverlapLevel_; level++) {
    OldGraph = OverlapGraph_;
    OldRowMap = OverlapRowMap_;

    OverlapImporter_ = (Epetra_Import *) OldGraph->Importer();
    OverlapRowMap_ = new Epetra_BlockMap(OverlapImporter_->TargetMap());

    if (level<OverlapLevel_)
      OverlapGraph_ = new Epetra_CrsGraph(Copy, *OverlapRowMap_, 0);
    else
      // On last iteration, we want to filter out all columns except those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlapGraph_ = new Epetra_CrsGraph(Copy, *OverlapRowMap_, *OverlapRowMap_, 0);

    EPETRA_CHK_ERR(OverlapGraph_->Import( *UserMatrixGraph, *OverlapImporter_, Insert));
    if (level<OverlapLevel_) {
      EPETRA_CHK_ERR(OverlapGraph_->TransformToLocal(DomainMap, RangeMap));
    }
    else {
      // Copy last OverlapImporter because we will use it later
      OverlapImporter_ = new Epetra_Import(*OverlapRowMap_, *DomainMap);
      EPETRA_CHK_ERR(OverlapGraph_->TransformToLocal(DomainMap, RangeMap));
    }

    if (level>1) {
      delete OldGraph;
      delete OldRowMap;
    }
  }

  return(0);
}

