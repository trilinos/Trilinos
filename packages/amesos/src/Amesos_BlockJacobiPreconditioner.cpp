#include "Amesos_ConfigDefs.h"
#include "Amesos_BlockJacobiPreconditioner.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"
#include <vector>
#include "Amesos_Partitioner.h"
#include "Amesos_Container.h"
#include "Amesos_ContainerFactory.h"
#include "Amesos_LocalRowMatrix.h"
#include "Amesos_InverseFactory.h"
#include "Amesos_Graph.h"
#include "Amesos_GraphRowMatrix.h"
#include "Amesos_Utils.h"

//==============================================================================
// FIXME: zero out all pointers
Amesos_BlockJacobiPreconditioner::
Amesos_BlockJacobiPreconditioner(char* Type,
					Epetra_RowMatrix* RowMatrix,
					Amesos_InverseFactory& Factory,
					int NumDomains,
					Teuchos::ParameterList& List) :
  Type_(Type),
  Factory_(Factory),
  RowMatrix_(RowMatrix),
  List_(List),
  NumLocalDomains_(NumDomains),
  IsPreconditionerComputed_(false),
  Graph_(0)
{

  // 0.- startup
 
  Label_ = string(Type) + " prec, " +
    "# local domains = " + Amesos_toString(NumLocalDomains());

  SubmatricesType_ = List_.get("submatrices type", "EpetraCrs");

  // 1.- sanity checks

  if (RowMatrix_->NumGlobalRows() != RowMatrix_->NumGlobalCols())
    AMESOS_CHK_ERRV(-1); // only square matrices

  // 3.- build graphs

  Graph_ = new Amesos_GraphRowMatrix(RowMatrix_);

  LocalizedMatrix_ = new Amesos_LocalRowMatrix(RowMatrix_);

  // 2.- extract the submatrices, and store them in an array

  AMESOS_CHK_ERRV(ExtractSubmatrices());
  
  //  .- compute the inverse 
  for (int i = 0 ; i < NumLocalDomains() ; ++i) {

    Containers_[i]->ComputeInverse((char*)Type_.c_str(),
				   Factory_, List_);

  }

  // 5.- that's all folks

  IsPreconditionerComputed_ = true;

}

//==============================================================================
Amesos_BlockJacobiPreconditioner::~Amesos_BlockJacobiPreconditioner()
{

  if (Graph_)
    delete Graph_;

  for (int i = 0 ; i < NumLocalDomains() ; ++i)
    delete Containers_[i];

}

//==============================================================================
int Amesos_BlockJacobiPreconditioner::ExtractSubmatrices()
{

  // build the partitioner
 
  // FIXME:
  List_.set("local parts", NumLocalDomains());
  List_.set("overlap level", (int)0);

  Partitioner_ = new Amesos_Partitioner(Graph_,(Amesos_Graph*)0, List_);

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
      int LRID = RowMatrix_->RowMatrixRowMap().LID(GRID);
      // set the GID of local row 
      Containers_[i]->GID(j) = LRID;
      
      int NumEntries;
      int Length = LocalizedMatrix_->MaxNumEntries();
      vector<double> Values;
      Values.resize(Length);
      vector<int> Indices;
      Indices.resize(Length);

      // FIXME: global or local indices??
      int ierr = 
	LocalizedMatrix_->ExtractMyRowCopy(LRID, Length, NumEntries, 
					   &Values[0], &Indices[0]);
      AMESOS_CHK_ERR(ierr);

      for (int k = 0 ; k < NumEntries ; ++k) {

	int LCID = Indices[k];
	int GCID = RowMatrix_->RowMatrixRowMap().GID(LCID);

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
int Amesos_BlockJacobiPreconditioner::
ApplyInverse(const Epetra_MultiVector& X, 
	     Epetra_MultiVector& Y) const
{

  if (IsPreconditionerComputed() == false)
    AMESOS_CHK_ERR(-4); // need to compute the prec first

  if (X.NumVectors() != Y.NumVectors())
    AMESOS_CHK_ERR(-1); // not valid


  Epetra_MultiVector Xtmp(X);

  Y.PutScalar(0.0);

  // cycle over all local subdomains

  for (int i = 0 ; i < NumLocalDomains() ; ++i) {

    int LID, GID;

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->GID(j);
      for (int k = 0 ; k < X.NumVectors() ; ++k) {
	Containers_[i]->RHS(j,k) = Xtmp[k][LID];
      }
    }

    AMESOS_CHK_ERR(Containers_[i]->ApplyInverse());

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->GID(j);
      for (int k = 0 ; k < X.NumVectors() ; ++k) {
	Y[k][LID] = Y[k][LID] + Containers_[i]->LHS(j,k);
      }
    }
  }

  return(0); 
}

//==============================================================================
char* Amesos_BlockJacobiPreconditioner::Label() const
{
  return(const_cast<char*>(Label_.c_str()));
}

//==============================================================================
int Amesos_BlockJacobiPreconditioner::
Apply(const Epetra_MultiVector& X, 
      Epetra_MultiVector& Y) const
{
  AMESOS_RETURN(RowMatrix_->Apply(X,Y));
}

//==============================================================================
const Epetra_Comm& Amesos_BlockJacobiPreconditioner::
Comm() const
{
  return(RowMatrix_->Comm());
}

//==============================================================================
const Epetra_Map& Amesos_BlockJacobiPreconditioner::
OperatorDomainMap() const
{
  return(RowMatrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map& Amesos_BlockJacobiPreconditioner::
OperatorRangeMap() const
{
  return(RowMatrix_->OperatorRangeMap());
}

