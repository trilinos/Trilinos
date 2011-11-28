//@HEADER
// ************************************************************************
// 
//              Ifpack_SubdomainFilter
//              Author: Radu Popescu <radu.popescu@epfl.ch>
//              Copyright 2011 EPFL - CMCS
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
// 3. Neither the name of the copyright holder nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY EPFL - CMCS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL EPFL - CMCS OR THE CONTRIBUTORS
// BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
// OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
//@HEADER

#ifdef IFPACK_SUBCOMM_CODE

#include <vector>
#include <algorithm>

#include "Ifpack_ConfigDefs.h"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Ifpack_SubdomainFilter.h"
#include "Ifpack_OverlappingRowMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_BLAS_wrappers.h"

using namespace Teuchos;

//==============================================================================
Ifpack_SubdomainFilter::Ifpack_SubdomainFilter(const RefCountPtr<const Epetra_RowMatrix>& Matrix,int subdomainId) :
  Matrix_(Matrix),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalCols_(0),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0)
{
  sprintf(Label_,"%s","Ifpack_SubdomainFilter");

  ImportVector_ = 0;
  ExportVector_ = 0;

  ovA_ = dynamic_cast<const Ifpack_OverlappingRowMatrix*>(&*Matrix_);

  if (ovA_) {
    Acrs_=dynamic_cast<const Epetra_CrsMatrix*>(&ovA_->A());
    NumMyRowsA_ = ovA_->A().NumMyRows();
    NumMyRowsB_ = ovA_->B().NumMyRows();
  } else {
    NumMyRowsA_ = Matrix->NumMyRows();
    NumMyRowsB_ = 0;
  }
  
#ifdef HAVE_MPI
  const Epetra_MpiComm *pComm = dynamic_cast<const Epetra_MpiComm*>( &(Matrix->Comm()) );
  assert(pComm != NULL);
  MPI_Comm_split(pComm->Comm(), subdomainId, pComm->MyPID(), &subdomainMPIComm_);
  SubComm_ = rcp( new Epetra_MpiComm(subdomainMPIComm_) );
#else
  SubComm_ = rcp( new Epetra_SerialComm );
#endif

  NumMyRows_ = Matrix->NumMyRows();
  SubComm_->SumAll(&NumMyRows_,&NumGlobalRows_,1);

  // Get the row GID corresponding to the process
  const Epetra_Map &globRowMap = Matrix->RowMatrixRowMap();
  std::vector<int> myRowsGID(NumMyRows_);
  globRowMap.MyGlobalElements(&myRowsGID[0]);

  // Gather the GID of all rows in the subdomain
  // Get the maximum number of local rows on each processor in the subdomain
  // Allocate a vector large enough (Proc cout * max number local rows
  // Gather all the local row indices. If a process has less than max num
  // local rows, pad the local vector with -1;
  // After the gather, eliminate the -1 elements from the target vector
  //
  // TODO: this section should be rewritten to avoid the GatherAll call.
  //       Right now it's not memory efficient. It could be moved to
  //       Epetra_Map since it a map operation

  int MaxLocalRows = 0;
  SubComm_->MaxAll(&NumMyRows_, &MaxLocalRows, 1);
  
  std::vector<int> SubdomainRowGID(SubComm_->NumProc() * MaxLocalRows);
  myRowsGID.resize(MaxLocalRows, -1);

  SubComm_->GatherAll(&myRowsGID[0],
                      &SubdomainRowGID[0],
                      MaxLocalRows);

  std::sort(SubdomainRowGID.begin(), SubdomainRowGID.end());

  int PaddingSize = SubdomainRowGID.size() - NumGlobalRows_;
  SubdomainRowGID.erase(SubdomainRowGID.begin(),
                         SubdomainRowGID.begin() + PaddingSize);

  // Get the col GID corresponding to the process
  const Epetra_Map &globColMap = Matrix->RowMatrixColMap();
  NumMyCols_ = globColMap.NumMyElements();
  std::vector<int> myGlobalCols(NumMyCols_);
  globColMap.MyGlobalElements(&myGlobalCols[0]);

  // Eliminate column GID that are outside the subdomain
  std::vector<int> filteredColGID;
  for (int i = 0; i < NumMyCols_; ++i) {
    if (std::find(SubdomainRowGID.begin(), SubdomainRowGID.end(),
      myGlobalCols[i]) != SubdomainRowGID.end()) {
      filteredColGID.push_back(myGlobalCols[i]);
    }
  }
  NumMyCols_ = filteredColGID.size();

  // Create the maps with the reindexed GIDs
  Map_ = rcp( new Epetra_Map(-1, NumMyRows_, &myRowsGID[0], globRowMap.IndexBase(), *SubComm_) );
  colMap_ = rcp( new Epetra_Map(-1, NumMyCols_, &filteredColGID[0], globColMap.IndexBase(), *SubComm_) );
  NumGlobalCols_ = NumGlobalRows_;

  NumEntries_.resize(NumMyRows_);
  MaxNumEntriesA_ = Matrix->MaxNumEntries();
  MaxNumEntries_ = Matrix->MaxNumEntries();

  Diagonal_ = rcp( new Epetra_Vector(*Map_) );
  if (Diagonal_ == Teuchos::null) IFPACK_CHK_ERRV(-5);

  Indices_.resize(MaxNumEntries_);
  Values_.resize(MaxNumEntries_);

  Ac_LIDMap_ = 0;
  Bc_LIDMap_ = 0;
  Ar_LIDMap_ = 0;
  Br_LIDMap_ = 0;

  if(ovA_){
    Ac_LIDMap_ = new int[ovA_->A().NumMyCols() + 1];
    for(int i = 0; i < ovA_->A().NumMyCols(); i++) {
      Ac_LIDMap_[i] = colMap_->LID(ovA_->A().RowMatrixColMap().GID(i));    
    }
    Bc_LIDMap_ = new int[ovA_->B().NumMyCols() + 1];
    for(int i = 0; i < ovA_->B().NumMyCols(); i++) {
      Bc_LIDMap_[i] = colMap_->LID(ovA_->B().RowMatrixColMap().GID(i));
    }
    Ar_LIDMap_ = new int[ovA_->A().NumMyRows() + 1];
    for(int i = 0; i < ovA_->A().NumMyRows(); i++) {
      Ar_LIDMap_[i] = Map_->LID(ovA_->A().RowMatrixRowMap().GID(i));    
    }
    Br_LIDMap_ = new int[ovA_->B().NumMyRows() + 1];
    for(int i = 0; i < ovA_->B().NumMyRows(); i++) {
      Br_LIDMap_[i] = Map_->LID(ovA_->B().RowMatrixRowMap().GID(i));
    }
  }

  int ActualMaxNumEntries = 0;

  for (int i = 0 ; i < NumMyRows_; ++i) {
    NumEntries_[i] = 0;
    int Nnz;
    IFPACK_CHK_ERRV(ExtractMyRowCopy(i, MaxNumEntries_, Nnz, &Values_[0], &Indices_[0]));

    for (int j = 0 ; j < Nnz; ++j) {
      if (Indices_[j] == i) {
        (*Diagonal_)[i] = Values_[j];
      }
    }

    if (Nnz > ActualMaxNumEntries) {
      ActualMaxNumEntries = Nnz;
    }

    NumMyNonzeros_ += Nnz;
    NumEntries_[i] = Nnz;
  }

  SubComm_->SumAll(&NumMyNonzeros_, &NumGlobalNonzeros_, 1);
  MaxNumEntries_ = ActualMaxNumEntries;

  int gpid = Matrix->Comm().MyPID();
  Exporter_ = null;
  Importer_ = null;

  if (!(RowMatrixRowMap().SameAs(OperatorRangeMap()))) {
    try{
      Exporter_ = rcp(new Epetra_Export(RowMatrixRowMap(), OperatorRangeMap()));
    }
    catch(...) {
      printf("** * gpid %d: Ifpack_SubdomainFilter ctor: problem creating Exporter_ * **\n\n",gpid);
    }
  }

  if (!(*colMap_).SameAs(*Map_)) {
    try{
      Importer_ = rcp(new Epetra_Import(*colMap_, *Map_));
    }
    catch(...) {
      printf("** * gpid %d: Ifpack_SubdomainFilter ctor: problem creating Importer_ * **\n\n",gpid);
    }
  }

} //Ifpack_SubdomainFilter() ctor

//==============================================================================
Ifpack_SubdomainFilter::~Ifpack_SubdomainFilter()
{
    if(Ac_LIDMap_) delete [] Ac_LIDMap_;
    if(Bc_LIDMap_) delete [] Bc_LIDMap_;
    if(Ar_LIDMap_) delete [] Ar_LIDMap_;
    if(Br_LIDMap_) delete [] Br_LIDMap_;

    if(ImportVector_) delete ImportVector_; 
} //Ifpack_SubdomainFilter destructor

//==============================================================================
int Ifpack_SubdomainFilter:: ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
    double *Values, int * Indices) const
{
  if ((MyRow < 0) || (MyRow >= NumMyRows_)) {
    IFPACK_CHK_ERR(-1);
  }

  assert(Length >= NumEntries_[MyRow]);

  int ierr;
  if (ovA_) {
    int LocRow = ovA_->A().RowMatrixRowMap().LID(Map_->GID(MyRow));      
    if (LocRow != -1) { 
      ierr = ovA_->A().ExtractMyRowCopy(LocRow, Length, NumEntries, Values, Indices);
      for(int j = 0;j < NumEntries; j++) {
        Indices[j] = Ac_LIDMap_[Indices[j]];
      }
    }
    else {
      LocRow = ovA_->B().RowMatrixRowMap().LID(Map_->GID(MyRow));      
      ierr = ovA_->B().ExtractMyRowCopy(LocRow, Length, NumEntries, Values, Indices);
      for(int j = 0; j < NumEntries; j++) {
        Indices[j] = Bc_LIDMap_[Indices[j]];
      }
    }
  } else {
    int Nnz = 0;
    ierr = Matrix_->ExtractMyRowCopy(MyRow, MaxNumEntriesA_, Nnz, &Values_[0], &Indices_[0]);
    IFPACK_CHK_ERR(ierr);

    NumEntries = 0;
    for (int j = 0 ; j < Nnz ; ++j) {
      int subdomainLocalIndex = colMap_->LID(Matrix_->RowMatrixColMap().GID(Indices_[j]));
      if (subdomainLocalIndex != -1) {
        Indices[NumEntries] = subdomainLocalIndex;
        Values[NumEntries] = Values_[j];
        NumEntries++;
      }
    }
  }

  return(ierr);
}

//==============================================================================
int Ifpack_SubdomainFilter::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  if (!Diagonal.Map().SameAs(*Map_)) {
    IFPACK_CHK_ERR(-1);
  }

  Diagonal = *Diagonal_;
  return(0);
}

//==============================================================================
int Ifpack_SubdomainFilter::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
    std::cout << "Ifpack_SubdomainFilter::Apply not implemented.\n";
    IFPACK_CHK_ERR(-1); // not implemented
}

//==============================================================================
void Ifpack_SubdomainFilter::UpdateImportVector(int NumVectors) const
{
}

//==============================================================================
void Ifpack_SubdomainFilter::UpdateExportVector(int NumVectors) const
{
}

//==============================================================================
int Ifpack_SubdomainFilter::ApplyInverse(const Epetra_MultiVector& X,
		 Epetra_MultiVector& Y) const
{
  std::cout << "Ifpack_SubdomainFilter::ApplyInverse not implemented.\n";
  IFPACK_CHK_ERR(-1); // not implemented
}

//==============================================================================
const Epetra_BlockMap& Ifpack_SubdomainFilter::Map() const
{
  return(*Map_);
}

#endif //ifdef IFPACK_SUBCOMM_CODE
