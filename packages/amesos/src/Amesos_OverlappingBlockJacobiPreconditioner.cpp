#include "Amesos_ConfigDefs.h"
#include "Amesos_OverlappingBlockJacobiPreconditioner.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Teuchos_ParameterList.hpp"
#include <vector>
#include "Amesos_Utils.h"
#include "Amesos_Partitioner.h"
#include "Amesos_Container.h"
#include "Amesos_ContainerFactory.h"
#include "Amesos_LocalRowMatrix.h"
#include "Amesos_InverseFactory.h"
#include "Amesos_Graph.h"
#include "Amesos_GraphEpetraCrs.h"

//==============================================================================
// FIXME: zero out all pointers
Amesos_OverlappingBlockJacobiPreconditioner::
Amesos_OverlappingBlockJacobiPreconditioner(char* Type,
					Epetra_CrsMatrix* Matrix,
					Amesos_InverseFactory& Factory,
					int OverlappingLevel,
					bool IsSymmetric,
					int NumDomains,
					Teuchos::ParameterList& List) :
  Type_(Type),
  Factory_(Factory),
  Matrix_(Matrix),
  List_(List),
  OverlappingLevel_(OverlappingLevel),
  NumLocalDomains_(NumDomains),
  IsSymmetric_(IsSymmetric_),
  IsPreconditionerComputed_(false),
  OverlappingMap_(0),
  OverlappingX_(0),
  OverlappingY_(0),
  Importer_(0),
  Exporter_(0),
  Graph_(0),
  OverlappingGraph_(0)
{

  // 0.- startup
 
  Label_ = string(Type) + " prec, ov = " +
    Amesos_toString(OverlappingLevel) +
    ", # local domains = " + Amesos_toString(NumLocalDomains());

  SubmatricesType_ = List_.get("submatrices type", "EpetraCrs");

  // 1.- sanity checks

  if (Matrix_->NumGlobalRows() != Matrix_->NumGlobalCols())
    AMESOS_CHK_ERRV(-1); // only square matrices

  // 2.- build the overlapping matrix

  if ((OverlappingLevel_ > 0) && (Comm().NumProc() > 1)) {

    OverlappingMatrix_ = Amesos_CreateOverlappingCrsMatrix(Matrix,
                                                           OverlappingLevel_);
    assert(OverlappingMatrix_ != 0);

  }
  else
    OverlappingMatrix_ = Matrix_;

  // 3.- build graphs

  Graph_ = new Amesos_GraphEpetraCrs(&Matrix_->Graph());
  if ((OverlappingLevel_ > 0) && (Comm().NumProc() > 1)) 
    OverlappingGraph_ = new Amesos_GraphEpetraCrs(&OverlappingMatrix_->Graph());
  else
    OverlappingGraph_ = Graph_;

  OverlappingMap_ = new Epetra_Map(OverlappingMatrix_->RowMatrixRowMap());

  // FIXME: a che servo??
  LocalizedOverlappingMatrix_ = new Amesos_LocalRowMatrix(OverlappingMatrix_);

  // 2.- extract the submatrices, and store them in an array

  AMESOS_CHK_ERRV(ExtractSubmatrices());
  
  // 3.- if SubmatricesType_ == "Epetra", use Amesos to build the
  //     local factorizations

  for (int i = 0 ; i < NumLocalDomains() ; ++i) {

    Containers_[i]->ComputeInverse((char*)Type_.c_str(),
				   Factory_, List_);

  }

  // 5.- free memory no longer needed

  /* FIXME
  if ((OverlappingLevel_ > 0) && (Comm().NumProc() > 1)) {
    delete OverlappingMatrix_;
    OverlappingMatrix_ = 0;
  }
  */

  Importer_ = new Epetra_Import(OverlappingMatrix_->RowMatrixRowMap(),
				Matrix->RowMatrixRowMap());
  Exporter_ = new Epetra_Export(Matrix->RowMatrixRowMap(),
				OverlappingMatrix_->RowMatrixRowMap());

  OverlappingX_ =
    new Epetra_MultiVector(OverlappingMatrix_->RowMatrixRowMap(),1);
  OverlappingY_ =
    new Epetra_MultiVector(OverlappingMatrix_->RowMatrixRowMap(),1);

  delete LocalizedOverlappingMatrix_;

  // 5.- that's all folks

  IsPreconditionerComputed_ = true;

}

//==============================================================================
Amesos_OverlappingBlockJacobiPreconditioner::
~Amesos_OverlappingBlockJacobiPreconditioner()
{

  if (Importer_)
    delete Importer_;

  if (Exporter_)
    delete Exporter_;

  if (OverlappingX_)
    delete OverlappingX_;

  if (OverlappingY_)
    delete OverlappingY_;

  if (OverlappingMap_)
    delete OverlappingMap_;

  if (Graph_)
    delete Graph_;

  for (int i = 0 ; i < NumLocalDomains() ; ++i)
    delete Containers_[i];

}

//==============================================================================
int Amesos_OverlappingBlockJacobiPreconditioner::ExtractSubmatrices()
{

  // build the partitioner
 
  // FIXME:
  List_.set("local parts", NumLocalDomains());
  List_.set("overlap level", OverlappingLevel());

  Partitioner_ = new Amesos_Partitioner(Graph_, OverlappingGraph_,
					List_);

  Partitioner_->Compute();

  NumLocalDomains_ = Partitioner_->NumLocalParts();

  Containers_.resize(NumLocalDomains());

  Amesos_ContainerFactory ContainerFactory;

  for (int i = 0 ; i < NumLocalDomains() ; ++i) {

    Containers_[i] = ContainerFactory.Create(SubmatricesType_);
    
    if (Containers_[i] == 0)
      AMESOS_CHK_ERR(-10);
    
    int rows = Partitioner_->NumRowsInPart(i);

    Containers_[i]->Shape(rows);

    // extract all rows corresponding to subgraph `i', and 
    // add the contributions for the subgraph

    for (int j = 0 ; j < rows ; ++j) {

      // Partitioner is returning the local ID
      int GRID = (*Partitioner_)(i,j);
      int LRID = OverlappingMatrix_->Graph().LRID(GRID);
      // set the GID of local row 
      Containers_[i]->GID(j) = LRID;
      
      int NumEntries;
      int Length = LocalizedOverlappingMatrix_->MaxNumEntries();
      vector<double> Values;
      Values.resize(Length);
      vector<int> Indices;
      Indices.resize(Length);

      // FIXME: global or local indices??
      int ierr = 
	LocalizedOverlappingMatrix_->ExtractMyRowCopy(LRID, Length, NumEntries, 
						      &Values[0], &Indices[0]);
      AMESOS_CHK_ERR(ierr);

      for (int k = 0 ; k < NumEntries ; ++k) {

	int LCID = Indices[k];
	int GCID = OverlappingMatrix_->Graph().GRID(LCID);

	int jj = -1;
	for (int kk = 0 ; kk < rows ; ++kk)
	  if ((*Partitioner_)(i,kk) == GCID)
	    jj = kk;

	if (jj != -1)
	  Containers_[i]->SetMatrixElement(j,jj,Values[k]);
	  
      }
    }

    Containers_[i]->Compute();

  }

  return(0);
}

//==============================================================================
int Amesos_OverlappingBlockJacobiPreconditioner::
ApplyInverse(const Epetra_MultiVector& X, 
	     Epetra_MultiVector& Y) const
{

  if (IsPreconditionerComputed() == false)
    AMESOS_CHK_ERR(-4); // need to compute the prec first

  if (X.NumVectors() != Y.NumVectors())
    AMESOS_CHK_ERR(-1); // not valid


  if (X.NumVectors() != OverlappingX_->NumVectors()) {
    delete OverlappingX_;
    delete OverlappingY_;
    OverlappingX_ =
      new Epetra_MultiVector(*OverlappingMap_, X.NumVectors());
    OverlappingY_ =
      new Epetra_MultiVector(*OverlappingMap_, X.NumVectors());
    assert (OverlappingX_ != 0);
    assert (OverlappingY_ != 0);
  }

  AMESOS_CHK_ERR(OverlappingX_->Export(X,*Exporter_,Insert));

  OverlappingY_->PutScalar(0.0);

  // cycle over all local subdomains

  for (int i = 0 ; i < NumLocalDomains() ; ++i) {

    int LID, GID;

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->GID(j);
      for (int k = 0 ; k < X.NumVectors() ; ++k) {
	Containers_[i]->RHS(j,k) = (*OverlappingX_)[k][LID];
      }
    }

    AMESOS_CHK_ERR(Containers_[i]->ApplyInverse());

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->GID(j);
      for (int k = 0 ; k < X.NumVectors() ; ++k) {
	(*OverlappingY_)[k][LID] = 
	  (*OverlappingY_)[k][LID] + Containers_[i]->LHS(j,k);
      }
    }
  }

  if (IsSymmetric()) {
    AMESOS_CHK_ERR(Y.Export(*OverlappingY_,*Importer_,Add));

  }
  else {
    AMESOS_CHK_ERR(Y.Export(*OverlappingY_,*Importer_,Zero));
  }

  return(0); 
}

//==============================================================================
char* Amesos_OverlappingBlockJacobiPreconditioner::Label() const
{
  return(const_cast<char*>(Label_.c_str()));
}

//==============================================================================
int Amesos_OverlappingBlockJacobiPreconditioner::
Apply(const Epetra_MultiVector& X, 
      Epetra_MultiVector& Y) const
{
  AMESOS_RETURN(Matrix_->Apply(X,Y));
}

//==============================================================================
const Epetra_Comm& Amesos_OverlappingBlockJacobiPreconditioner::
Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map& Amesos_OverlappingBlockJacobiPreconditioner::
OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map& Amesos_OverlappingBlockJacobiPreconditioner::
OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}
