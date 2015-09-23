/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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

