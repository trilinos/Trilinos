/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * October 20, 2002, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "Ifpack_OverlapGraph.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"

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

    OverlapGraph_ = new Epetra_CrsGraph(Copy, *OverlapRowMap_, 0);
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

