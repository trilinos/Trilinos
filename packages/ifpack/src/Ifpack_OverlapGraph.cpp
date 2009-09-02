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

#include <Teuchos_ParameterList.hpp>
#include <ifp_parameters.h>

//==============================================================================
Ifpack_OverlapGraph::Ifpack_OverlapGraph(const Teuchos::RefCountPtr<const Epetra_CrsGraph>& UserMatrixGraph_in, int OverlapLevel_in)
  : UserMatrixGraph_(UserMatrixGraph_in),
    OverlapLevel_(OverlapLevel_in)
{
  // Test for non-trivial overlap here so we can use it later.
  IsOverlapped_ = (OverlapLevel_in>0 && UserMatrixGraph_in->DomainMap().DistributedGlobal());

  ConstructOverlapGraph(UserMatrixGraph_in);

}
//==============================================================================
Ifpack_OverlapGraph::Ifpack_OverlapGraph(const Teuchos::RefCountPtr<const Epetra_RowMatrix>& UserMatrix_in, int OverlapLevel_in)
  : UserMatrix_(UserMatrix_in),
    OverlapLevel_(OverlapLevel_in)
{
  // Test for non-trivial overlap here so we can use it later.
  IsOverlapped_ = (OverlapLevel_in>0 && UserMatrix_in->OperatorDomainMap().DistributedGlobal());
  
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
    if (OverlapGraph_!=Teuchos::null) OverlapGraph_ = Teuchos::rcp( new Epetra_CrsGraph(*OverlapGraph_) );
    if (OverlapRowMap_!=Teuchos::null) OverlapRowMap_ = Teuchos::rcp( new Epetra_BlockMap(*OverlapRowMap_) );
  }
  
}

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

//==============================================================================
int Ifpack_OverlapGraph::ConstructOverlapGraph(const Teuchos::RefCountPtr<const Epetra_CrsGraph>& UserMatrixGraph) {

  if (!IsOverlapped_) {
    OverlapGraph_ = Teuchos::rcp_const_cast<Epetra_CrsGraph>( UserMatrixGraph );
    OverlapRowMap_ = Teuchos::rcp( (Epetra_BlockMap *) &UserMatrixGraph->RowMap(), false );
    return(0); // All done
  }

  Teuchos::RefCountPtr<Epetra_CrsGraph> OldGraph;
  Teuchos::RefCountPtr<Epetra_BlockMap> OldRowMap;
  const Epetra_BlockMap DomainMap = UserMatrixGraph->DomainMap();
  const Epetra_BlockMap RangeMap = UserMatrixGraph->RangeMap();

  for (int level=1; level <= OverlapLevel_; level++) {
    OldGraph = OverlapGraph_;
    OldRowMap = OverlapRowMap_;

    OverlapImporter_ = Teuchos::rcp( (Epetra_Import *) OldGraph->Importer(), false );
    OverlapRowMap_ = Teuchos::rcp( new Epetra_BlockMap(OverlapImporter_->TargetMap()) );

    if (level<OverlapLevel_)
      OverlapGraph_ = Teuchos::rcp( new Epetra_CrsGraph(Copy, *OverlapRowMap_, 0) );
    else
      // On last iteration, we want to filter out all columns except those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlapGraph_ = Teuchos::rcp( new Epetra_CrsGraph(Copy, *OverlapRowMap_, *OverlapRowMap_, 0) );

    EPETRA_CHK_ERR(OverlapGraph_->Import( *UserMatrixGraph, *OverlapImporter_, Insert));
    if (level<OverlapLevel_) {
      EPETRA_CHK_ERR(OverlapGraph_->FillComplete(DomainMap, RangeMap));
    }
    else {
      // Copy last OverlapImporter because we will use it later
      OverlapImporter_ = Teuchos::rcp( new Epetra_Import(*OverlapRowMap_, DomainMap) );
      EPETRA_CHK_ERR(OverlapGraph_->FillComplete(DomainMap, RangeMap));
    }

  }

  return(0);
}

