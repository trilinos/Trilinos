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

#ifdef IFPACK_NODE_AWARE_CODE

#include "Ifpack_ConfigDefs.h"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Ifpack_NodeFilter.h"
#include "Ifpack_OverlappingRowMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"

using namespace Teuchos;

extern int ML_NODE_ID;

//==============================================================================
Ifpack_NodeFilter::Ifpack_NodeFilter(const RefCountPtr<const Epetra_RowMatrix>& Matrix,int nodeID) :
  Matrix_(Matrix),
  NumMyRows_(0),
  NumGlobalRows_(0),
  NumMyNonzeros_(0),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0)
{
  sprintf(Label_,"%s","Ifpack_NodeFilter");

  printf("------ --- - Entering Ifpack_NodeFilter::Ifpack_NodeFilter ctor\n"); fflush(stdout);
  ImportVector_=null;
  //ExportVector_=null;
  ExportVector_=0;

#ifdef HAVE_MPI
  const Epetra_MpiComm *pComm = dynamic_cast<const Epetra_MpiComm*>( &(Matrix->Comm()) );
  assert(pComm != NULL);
  MPI_Comm_split(pComm->Comm(),nodeID,pComm->MyPID(),&nodeMPIComm_);
  SubComm_ = rcp( new Epetra_MpiComm(nodeMPIComm_) );
#else
  SubComm_ = rcp( new Epetra_SerialComm );
#endif

  NumMyRows_ = Matrix->NumMyRows();
  SubComm_->SumAll(&NumMyRows_,&NumGlobalRows_,1);
  NumMyCols_ = Matrix->NumMyCols();

  // build row map, based on the local communicator
  try {
    const Epetra_Map &globRowMap = Matrix->RowMatrixRowMap();
    int *myGlobalElts =  globRowMap.MyGlobalElements();
    int numMyElts = globRowMap.NumMyElements();
    Map_ = rcp( new Epetra_Map(-1,numMyElts,myGlobalElts,globRowMap.IndexBase(),*SubComm_) );
  }
  catch(...) {
    printf("** * Ifpack_NodeFilter ctor: problem creating row map * **\n\n");
  }


  // build the column map, but don't use a copy constructor, b/c local communicator SubComm_ is
  // different from that of Matrix.
  try {
    const Epetra_Map &globColMap = Matrix->RowMatrixColMap();
    int *myGlobalElts =  globColMap.MyGlobalElements();
    int numMyElts = globColMap.NumMyElements();
    colMap_ = rcp( new Epetra_Map(-1,numMyElts,myGlobalElts,globColMap.IndexBase(),*SubComm_) );
  }
  catch(...) {
    printf("** * Ifpack_NodeFilter ctor: problem creating col map * **\n\n");
  }

  // NumEntries_ will contain the actual number of nonzeros
  // for each localized row (that is, without external nodes,
  // and always with the diagonal entry)
  NumEntries_.resize(NumMyRows_);

  // want to store the diagonal vector. FIXME: am I really useful?
  Diagonal_ = rcp( new Epetra_Vector(*Map_) );
  if (Diagonal_ == Teuchos::null) IFPACK_CHK_ERRV(-5);

  // store this for future access to ExtractMyRowCopy().
  // This is the # of nonzeros in the non-local matrix
  MaxNumEntriesA_ = Matrix->MaxNumEntries();
  // tentative value for MaxNumEntries. This is the number of
  // nonzeros in the local matrix
  MaxNumEntries_ = Matrix->MaxNumEntries();

  // ExtractMyRowCopy() will use these vectors
  Indices_.resize(MaxNumEntries_);
  Values_.resize(MaxNumEntries_);

  // now compute:
  // - the number of nonzero per row
  // - the total number of nonzeros
  // - the diagonal entries

  // compute nonzeros (total and per-row), and store the
  // diagonal entries (already modified)
  int ActualMaxNumEntries = 0;

  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    
    NumEntries_[i] = 0;
    int Nnz, NewNnz = 0;
    IFPACK_CHK_ERRV(ExtractMyRowCopy(i,MaxNumEntries_,Nnz,&Values_[0],&Indices_[0]));

    for (int j = 0 ; j < Nnz ; ++j) {
      NewNnz++;
      if (Indices_[j] == i) (*Diagonal_)[i] = Values_[j];
    }

    if (NewNnz > ActualMaxNumEntries)
      ActualMaxNumEntries = NewNnz;

    NumMyNonzeros_ += NewNnz;
    NumEntries_[i] = NewNnz;

  }

  SubComm_->SumAll(&NumMyNonzeros_,&NumGlobalNonzeros_,1);
  MaxNumEntries_ = ActualMaxNumEntries;

  int gpid = Matrix->Comm().MyPID();
  int lpid = SubComm_->MyPID();

  Exporter_ = null;
  Importer_ = null;
  // Check if non-trivial import/export operators
  if (!(RowMatrixRowMap().SameAs(OperatorRangeMap()))) {
    try{Exporter_ = rcp(new Epetra_Export(RowMatrixRowMap(), OperatorRangeMap()));}
    catch(...) {
      printf("** * gpid %d: Ifpack_NodeFilter ctor: problem creating Exporter_ * **\n\n",gpid);
    }
  }
  //if (gpid > 1) sleep(8);
/*
  if (gpid == 0)
    printf(">>> node 0 <<<\n");
  if (gpid == 2)
    printf(">>> node 1 <<<\n");
  if (lpid == 0)
    printf("=======================================\ntarget: RowMatrixColMap()\n====================================\n");
  cout << RowMatrixColMap() << endl;
  sleep(1);
  if (gpid == 0)
    printf("=======================================\nsource: OperatorDomainMap()\n=======================================\n");
  cout << OperatorDomainMap() << endl;
  sleep(1);
*/
/*
  if (!(RowMatrixColMap().SameAs(OperatorDomainMap()))) {
    //TODO change this to RCP
    try{Importer_ = new Epetra_Import(RowMatrixColMap(), OperatorDomainMap());}
    catch(...) {
      printf("** * gpid %d: Ifpack_NodeFilter ctor: problem creating Importer_ * **\n\n",gpid);
    }
*/
  if (!(*colMap_).SameAs(*Map_)) {
    //TODO change this to RCP
    try{Importer_ = rcp(new Epetra_Import(*colMap_, *Map_));}
    catch(...) {
      printf("** * gpid %d: Ifpack_NodeFilter ctor: problem creating Importer_ * **\n\n",gpid);
    }
/*
    if (lpid == 0)
    printf("=======================================\nIfpack_NodeFilter Importer_ on node %d\n=======================================\n",(gpid == 0 ? 0: 1)); fflush(stdout);
    cout << *Importer_ << endl;
    if (lpid == 0)
    printf("=======================================\nIfpack_NodeFilter colmap on node %d\n=======================================\n",(gpid == 0 ? 0: 1)); fflush(stdout);
    cout << *colMap_ << endl;
*/
  }

} //Ifpack_NodeFilter() ctor

//==============================================================================
int Ifpack_NodeFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  if ((MyRow < 0) || (MyRow >= NumMyRows_)) {
    IFPACK_CHK_ERR(-1); // range not valid
  }

  if (Length < NumEntries_[MyRow])
    assert(1==0);
    //return(-1);

  //TODO  will we have problems in the original use case?
  // always extract using the object Values_ and Indices_.
  // This is because I need more space than that given by
  // the user (for the external nodes)

  int GlobRow = Map_->GID(MyRow);
  //TODO this is really bad to do here from a performance standpoint
  const Ifpack_OverlappingRowMatrix *ovA = dynamic_cast<const Ifpack_OverlappingRowMatrix*>(&*Matrix_);
  //int ierr = ovA->ExtractGlobalRowCopy(GlobRow,MaxNumEntriesA_,&NumEntries, Values,Indices);
  int ierr = ovA->ExtractGlobalRowCopy(GlobRow,Length,NumEntries, Values,Indices);
  assert(ierr==0);

  // populate the user's vectors, add diagonal if not found

  for (int j = 0 ; j < NumEntries ; ++j) {
    // only local indices
      //printf("   MyRow = %d,   Indices_[%d] = %d, Values_[%d] = %g\n",MyRow, j, Indices_[j], j, Values_[j]);
      //fflush(stdout);
      Indices[j] = colMap_->LID(Indices[j]);
      assert(Indices[j] != -1);
  }
    
  return(ierr);

}

//==============================================================================
int Ifpack_NodeFilter::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  if (!Diagonal.Map().SameAs(*Map_))
    IFPACK_CHK_ERR(-1);
  Diagonal = *Diagonal_;
  return(0);
}

//==============================================================================
/*
//old Apply (no communication)
int Ifpack_NodeFilter::Apply(const Epetra_MultiVector& X,
	  Epetra_MultiVector& Y) const 
{

  // skip expensive checks, I suppose input data are ok

  Y.PutScalar(0.0);
  int NumVectors = Y.NumVectors();

  double** X_ptr;
  double** Y_ptr;
  X.ExtractView(&X_ptr);
  Y.ExtractView(&Y_ptr);

  for (int i = 0 ; i < NumRows_ ; ++i) {
    
    int Nnz;
    int ierr = Matrix_->ExtractMyRowCopy(i,MaxNumEntriesA_,Nnz,&Values_[0],
                                         &Indices_[0]);
    IFPACK_CHK_ERR(ierr);

    for (int j = 0 ; j < Nnz ; ++j) {
      if (Indices_[j] < NumRows_ ) {
	for (int k = 0 ; k < NumVectors ; ++k)
	  Y_ptr[k][i] += Values_[j] * X_ptr[k][Indices_[j]];
      }
    }
  }

  return(0);
}
*/
int Ifpack_NodeFilter::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  //
  // This function forms the product Y = A * X.
  //

  int NumEntries;

  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) {
    EPETRA_CHK_ERR(-1); // Need same number of vectors in each MV
  }

  // Make sure Import and Export Vectors are compatible
  UpdateImportVector(NumVectors);
  UpdateExportVector(NumVectors);

  double ** Xp = (double**) X.Pointers();
  double ** Yp = (double**) Y.Pointers();


  // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
  if (Importer()!=0) {
    EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
    Xp = (double**)ImportVector_->Pointers();
  }

  // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
  if (Exporter()!=0) {
    Yp = (double**)ExportVector_->Pointers();
  }

  // Do actual computation
  for(int i = 0; i < NumMyRows_; i++) {
    EPETRA_CHK_ERR(ExtractMyRowCopy(i, MaxNumEntries(), NumEntries, &Values_[0], &Indices_[0]));
    for (int k=0; k<NumVectors; k++) {
      double sum = 0.0;
      for(int j = 0; j < NumEntries; j++)
        sum += Values_[j]*Xp[k][Indices_[j]];
      Yp[k][i] = sum;
    }
  }
  
  if (Exporter()!=0) {
    Y.PutScalar(0.0);  // Make sure target is zero
    Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
  }
  // Handle case of rangemap being a local replicated map
  if (!OperatorRangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());

  return(0);
} //Apply

//==============================================================================

void Ifpack_NodeFilter::UpdateImportVector(int NumVectors) const {
  if(ImportVector_ == null || ImportVector_->NumVectors() != NumVectors)
    ImportVector_ = rcp(new Epetra_MultiVector(Importer_->TargetMap(),NumVectors));
}

//=======================================================================================================

void Ifpack_NodeFilter::UpdateExportVector(int NumVectors) const {
  if(Exporter() != 0) {
    if(ExportVector_ != 0) {
      if(ExportVector_->NumVectors() != NumVectors) {
     delete ExportVector_;
     ExportVector_= 0;
      }
    }
    if(ExportVector_ == 0)
      ExportVector_ = new Epetra_MultiVector(Exporter_->SourceMap(),NumVectors); // Create Export vector if needed
  }
  return;
/*
  if(ExportVector_ == null || ExportVector_->NumVectors() != NumVectors)
    ExportVector_ = rcp(new Epetra_MultiVector(Exporter_->SourceMap(),NumVectors));
*/
}

//=======================================================================================================
int Ifpack_NodeFilter::ApplyInverse(const Epetra_MultiVector& X,
		 Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-1); // not implemented
}

//==============================================================================
const Epetra_BlockMap& Ifpack_NodeFilter::Map() const
{
  return(*Map_);
}

#endif //ifdef IFPACK_NODE_AWARE_CODE
