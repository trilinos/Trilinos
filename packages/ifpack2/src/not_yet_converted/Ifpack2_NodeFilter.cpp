// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef IFPACK2_NODE_AWARE_CODE

#include "Ifpack2_ConfigDefs.hpp"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_IntVector.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_BlockMap.hpp"
#include "Ifpack2_NodeFilter.hpp"
#include "Ifpack2_OverlappingRowMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_CrsMatrix.hpp"

using namespace Teuchos;

extern int ML_NODE_ID;

//==============================================================================
Ifpack2_NodeFilter::Ifpack2_NodeFilter(const RefCountPtr<const Tpetra_RowMatrix>& Matrix,int nodeID) :
  Matrix_(Matrix),
  NumMyRows_(0),
  NumGlobalRows_(0),
  NumMyNonzeros_(0),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0)
{
  sprintf(Label_,"%s","Ifpack2_NodeFilter");

  //ImportVector_=null;
  //ExportVector_=null;
  ImportVector_=0;
  ExportVector_=0;

  ovA_ = dynamic_cast<const Ifpack2_OverlappingRowMatrix*>(&*Matrix_);
  assert(ovA_ != 0);

  Acrs_=dynamic_cast<const Tpetra_CrsMatrix*>(&ovA_->A());

  NumMyRowsA_ = ovA_->A().NumMyRows();
  NumMyRowsB_ = ovA_->B().NumMyRows();
  
#ifdef HAVE_MPI
  const Tpetra_MpiComm *pComm = dynamic_cast<const Tpetra_MpiComm*>( &(Matrix->Comm()) );
  assert(pComm != NULL);
  MPI_Comm_split(pComm->Comm(),nodeID,pComm->MyPID(),&nodeMPIComm_);
  SubComm_ = rcp( new Tpetra_MpiComm(nodeMPIComm_) );
#else
  SubComm_ = rcp( new Tpetra_SerialComm );
#endif

  NumMyRows_ = Matrix->NumMyRows();
  SubComm_->SumAll(&NumMyRows_,&NumGlobalRows_,1);
  NumMyCols_ = Matrix->NumMyCols();

  // build row map, based on the local communicator
  try {
    const Tpetra_Map &globRowMap = Matrix->RowMatrixRowMap();
    int *myGlobalElts =  globRowMap.MyGlobalElements();
    int numMyElts = globRowMap.NumMyElements();
    Map_ = rcp( new Tpetra_Map(-1,numMyElts,myGlobalElts,globRowMap.IndexBase(),*SubComm_) );
  }
  catch(...) {
    printf("** * Ifpack2_NodeFilter ctor: problem creating row map * **\n\n");
  }


  // build the column map, but don't use a copy constructor, b/c local communicator SubComm_ is
  // different from that of Matrix.
  try {
    const Tpetra_Map &globColMap = Matrix->RowMatrixColMap();
    int *myGlobalElts =  globColMap.MyGlobalElements();
    int numMyElts = globColMap.NumMyElements();
    colMap_ = rcp( new Tpetra_Map(-1,numMyElts,myGlobalElts,globColMap.IndexBase(),*SubComm_) );
  }
  catch(...) {
    printf("** * Ifpack2_NodeFilter ctor: problem creating col map * **\n\n");
  }

  // NumEntries_ will contain the actual number of nonzeros
  // for each localized row (that is, without external nodes,
  // and always with the diagonal entry)
  NumEntries_.resize(NumMyRows_);

  // want to store the diagonal vector. FIXME: am I really useful?
  Diagonal_ = rcp( new Tpetra_Vector(*Map_) );
  if (Diagonal_ == Teuchos::null) IFPACK2_CHK_ERRV(-5);

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


  // CMS: [A|B]-Local to Overlap-Local Column Indices
  if(ovA_){
    Ac_LIDMap_=new int[ovA_->A().NumMyCols()+1];
    for(int i=0;i<ovA_->A().NumMyCols();i++) Ac_LIDMap_[i]=colMap_->LID(ovA_->A().RowMatrixColMap().GID(i));    
    Bc_LIDMap_=new int[ovA_->B().NumMyCols()+1];
    for(int i=0;i<ovA_->B().NumMyCols();i++) Bc_LIDMap_[i]=colMap_->LID(ovA_->B().RowMatrixColMap().GID(i));

    Ar_LIDMap_=new int[ovA_->A().NumMyRows()+1];
    for(int i=0;i<ovA_->A().NumMyRows();i++) Ar_LIDMap_[i]=Map_->LID(ovA_->A().RowMatrixRowMap().GID(i));    
    Br_LIDMap_=new int[ovA_->B().NumMyRows()+1];
    for(int i=0;i<ovA_->B().NumMyRows();i++) Br_LIDMap_[i]=Map_->LID(ovA_->B().RowMatrixRowMap().GID(i));

  }
  // end CMS


  // compute nonzeros (total and per-row), and store the
  // diagonal entries (already modified)
  int ActualMaxNumEntries = 0;

  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    
    NumEntries_[i] = 0;
    int Nnz, NewNnz = 0;
    IFPACK2_CHK_ERRV(ExtractMyRowCopy(i,MaxNumEntries_,Nnz,&Values_[0],&Indices_[0]));

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
    try{Exporter_ = rcp(new Tpetra_Export(RowMatrixRowMap(), OperatorRangeMap()));}
    catch(...) {
      printf("** * gpid %d: Ifpack2_NodeFilter ctor: problem creating Exporter_ * **\n\n",gpid);
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
    try{Importer_ = new Tpetra_Import(RowMatrixColMap(), OperatorDomainMap());}
    catch(...) {
      printf("** * gpid %d: Ifpack2_NodeFilter ctor: problem creating Importer_ * **\n\n",gpid);
    }
*/



  if (!(*colMap_).SameAs(*Map_)) {
    //TODO change this to RCP
    try{Importer_ = rcp(new Tpetra_Import(*colMap_, *Map_));}
    catch(...) {
      printf("** * gpid %d: Ifpack2_NodeFilter ctor: problem creating Importer_ * **\n\n",gpid);
    }
/*
    if (lpid == 0)
    printf("=======================================\nIfpack2_NodeFilter Importer_ on node %d\n=======================================\n",(gpid == 0 ? 0: 1)); fflush(stdout);
    cout << *Importer_ << endl;
    if (lpid == 0)
    printf("=======================================\nIfpack2_NodeFilter colmap on node %d\n=======================================\n",(gpid == 0 ? 0: 1)); fflush(stdout);
    cout << *colMap_ << endl;
*/
  }

} //Ifpack2_NodeFilter() ctor

//==============================================================================
int Ifpack2_NodeFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  if ((MyRow < 0) || (MyRow >= NumMyRows_)) {
    IFPACK2_CHK_ERR(-1); // range not valid
  }

  if (Length < NumEntries_[MyRow])
    assert(1==0);
    //return(-1);

  //TODO  will we have problems in the original use case?
  // always extract using the object Values_ and Indices_.
  // This is because I need more space than that given by
  // the user (for the external nodes)

  int ierr;
  if (ovA_) {
    int LocRow=ovA_->A().RowMatrixRowMap().LID(Map_->GID(MyRow));      
    if (LocRow != -1) { 
      ierr=ovA_->A().ExtractMyRowCopy(LocRow,Length,NumEntries,Values,Indices);
      for(int j=0;j<NumEntries;j++)
	Indices[j]=Ac_LIDMap_[Indices[j]];
      
    }
    else {
      LocRow=ovA_->B().RowMatrixRowMap().LID(Map_->GID(MyRow));      
      ierr=ovA_->B().ExtractMyRowCopy(LocRow,Length,NumEntries,Values,Indices);
      for(int j=0;j<NumEntries;j++)
	Indices[j]=Bc_LIDMap_[Indices[j]];
    }
  } else {
    int Nnz;
    ierr = Matrix_->ExtractMyRowCopy(MyRow,MaxNumEntriesA_,Nnz, &Values_[0],&Indices_[0]);
    IFPACK2_CHK_ERR(ierr);

    NumEntries = 0;
    for (int j = 0 ; j < Nnz ; ++j) {
      // only local indices
      if (Indices_[j] < NumMyRows_ ) {
        Indices[NumEntries] = Indices_[j];
        Values[NumEntries] = Values_[j];
        ++NumEntries;
      }
    }
  }
    
  return(ierr);

}

//==============================================================================
int Ifpack2_NodeFilter::ExtractDiagonalCopy(Tpetra_Vector & Diagonal) const
{
  if (!Diagonal.Map().SameAs(*Map_))
    IFPACK2_CHK_ERR(-1);
  Diagonal = *Diagonal_;
  return(0);
}

//==============================================================================
/*
//old Apply (no communication)
int Ifpack2_NodeFilter::Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
	  Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const 
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
    IFPACK2_CHK_ERR(ierr);

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
int Ifpack2_NodeFilter::Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const {
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
  assert(ovA_!=0);
  int *MyRows;
  double *MyValues;
  int *MyIndices;

  if(Acrs_){
    // A rows - CrsMatrix Case
    IFPACK2_CHK_ERR(Acrs_->ExtractCrsDataPointers(MyRows,MyIndices,MyValues));
    //special case NumVectors==1
    if (NumVectors==1) {
      for(int i=0;i<NumMyRowsA_;i++) {
        int LocRow=Ar_LIDMap_[i];
        double sum = 0.0;
        for(int j = MyRows[i]; j < MyRows[i+1]; j++)
          sum += MyValues[j]*Xp[0][Ac_LIDMap_[MyIndices[j]]];          
        Yp[0][LocRow] = sum;
      }
    }
    else {
      for(int i=0;i<NumMyRowsA_;i++) {
        int LocRow=Ar_LIDMap_[i];
        for (int k=0; k<NumVectors; k++) {
          double sum = 0.0;
          for(int j = MyRows[i]; j < MyRows[i+1]; j++)
            sum += MyValues[j]*Xp[k][Ac_LIDMap_[MyIndices[j]]];          
          Yp[k][LocRow] = sum;
        }
      }
    }
  }
  else{
    // A rows - RowMatrix Case
    MyValues=&Values_[0];
    MyIndices=&Indices_[0];
    for(int i=0;i<NumMyRowsA_;i++) {
      ovA_->A().ExtractMyRowCopy(i,MaxNumEntries_,NumEntries,MyValues,MyIndices);
      int LocRow=Ar_LIDMap_[i];
      for (int k=0; k<NumVectors; k++) {
        double sum = 0.0;
        for(int j = 0; j < NumEntries; j++)
        sum += MyValues[j]*Xp[k][Ac_LIDMap_[MyIndices[j]]];          
        Yp[k][LocRow] = sum;
      }
    }
  }

  // B rows, always CrsMatrix
  IFPACK2_CHK_ERR(ovA_->B().ExtractCrsDataPointers(MyRows,MyIndices,MyValues));
  //special case NumVectors==1
  if (NumVectors==1) {
    for(int i=0;i<NumMyRowsB_;i++) {
      int LocRow=Br_LIDMap_[i];
      double sum = 0.0;
      for(int j = MyRows[i]; j < MyRows[i+1]; j++)
        sum += MyValues[j]*Xp[0][Bc_LIDMap_[MyIndices[j]]];          
      Yp[0][LocRow] = sum;
    }
  } else {
    for(int i=0;i<NumMyRowsB_;i++) {
      int LocRow=Br_LIDMap_[i];
      for (int k=0; k<NumVectors; k++) { //FIXME optimization, check for NumVectors=1
        double sum = 0.0;
        for(int j = MyRows[i]; j < MyRows[i+1]; j++)
          sum += MyValues[j]*Xp[k][Bc_LIDMap_[MyIndices[j]]];          
        Yp[k][LocRow] = sum;
      }
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

void Ifpack2_NodeFilter::UpdateImportVector(int NumVectors) const {
  if(Importer() != 0) {
    if(ImportVector_ != 0) {
      if(ImportVector_->NumVectors() != NumVectors) {
     delete ImportVector_;
     ImportVector_= 0;
      }
    }
    if(ImportVector_ == 0)
      ImportVector_ = new Tpetra_MultiVector(Importer_->TargetMap(),NumVectors); // Create Import vector if needed
  }
  return;
/*
  if(ImportVector_ == null || ImportVector_->NumVectors() != NumVectors)
    ImportVector_ = rcp(new Tpetra_MultiVector(Importer_->TargetMap(),NumVectors));
*/
}

//=======================================================================================================

void Ifpack2_NodeFilter::UpdateExportVector(int NumVectors) const {
  if(Exporter() != 0) {
    if(ExportVector_ != 0) {
      if(ExportVector_->NumVectors() != NumVectors) {
     delete ExportVector_;
     ExportVector_= 0;
      }
    }
    if(ExportVector_ == 0)
      ExportVector_ = new Tpetra_MultiVector(Exporter_->SourceMap(),NumVectors); // Create Export vector if needed
  }
  return;
/*
  if(ExportVector_ == null || ExportVector_->NumVectors() != NumVectors)
    ExportVector_ = rcp(new Tpetra_MultiVector(Exporter_->SourceMap(),NumVectors));
*/
}

//=======================================================================================================
int Ifpack2_NodeFilter::ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
		 Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  IFPACK2_CHK_ERR(-1); // not implemented
}

//==============================================================================
const Tpetra_BlockMap& Ifpack2_NodeFilter::Map() const
{
  return(*Map_);

}

#endif //ifdef IFPACK2_NODE_AWARE_CODE
