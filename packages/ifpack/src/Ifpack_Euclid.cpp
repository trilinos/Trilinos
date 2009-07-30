/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
#include "Ifpack_Euclid.h"
#if defined(HAVE_EUCLID) && defined(HAVE_MPI)

#include "Ifpack_Utils.h"
#include <algorithm>
#include "Epetra_MpiComm.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "getRow_dh.h"
#ifdef IFPACK_NODE_AWARE_CODE
extern int ML_NODE_ID;
#endif

using Teuchos::RCP;
using Teuchos::rcp;

Ifpack_Euclid::Ifpack_Euclid(Epetra_CrsMatrix* A):
  A_(rcp(A,false)),
  UseTranspose_(false),
  IsInitialized_(false),
  IsComputed_(false),
  Label_(),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  Time_(A_->Comm()),
  NiceRowMap_(true),
  SetLevel_(1),
  SetBJ_(0),
  SetStats_(0),
  SetMem_(0),
  SetSparse_(0.0),
  SetRowScale_(0),
  SetILUT_(0.0)
{
  // First make a copy of the Epetra Matrix as a Euclid Matrix
  int ilower = A_->RowMap().MinMyGID();
  int iupper = A_->RowMap().MaxMyGID();
  // Need to check if the RowMap is the way Euclid expects (if not more difficult)
  std::vector<int> ilowers; ilowers.resize(Comm().NumProc());
  std::vector<int> iuppers; iuppers.resize(Comm().NumProc());
  int myLower[1]; myLower[0] = ilower;
  int myUpper[1]; myUpper[0] = iupper;
  Comm().GatherAll(myLower, &ilowers[0], 1);
  Comm().GatherAll(myUpper, &iuppers[0], 1);
  for(int i = 0; i < Comm().NumProc()-1; i++){
    NiceRowMap_ = (NiceRowMap_ && iuppers[i]+1 == ilowers[i+1]);
  }
  if(!NiceRowMap_){
    ilower = (A_->NumGlobalRows() / Comm().NumProc())*Comm().MyPID();
    iupper = (A_->NumGlobalRows() / Comm().NumProc())*(Comm().MyPID()+1)-1;
    if(Comm().MyPID() == Comm().NumProc()-1){
      iupper = A_-> NumGlobalRows()-1;
    }
  }

  std::vector<int> rows; rows.resize(iupper - ilower +1);
  for(int i = ilower; i <= iupper; i++){
    rows[i-ilower] = i;
  }
  MySimpleMap_ = rcp(new Epetra_Map(-1, iupper-ilower+1, &rows[0], 0, Comm()));
  for(int i = 0; i < A_->NumMyRows(); i++){
    int *indices;
    int len;
    A_->Graph().ExtractMyRowView(i, len, indices);
    for(int j = 0; j < len; j++){
      indices[j] = A_->GCID(indices[j]);
    }
  }
}

void Ifpack_Euclid::Destroy(){
  if(IsComputed()){
    Euclid_dhDestroy(eu);
  }
  if(IsInitialized()){
    Mem_dhDestroy(mem_dh);
    mem_dh = NULL;
  }
  for(int i = 0; i < A_->NumMyRows(); i++){
    int *indices;
    int len;
    A_->Graph().ExtractMyRowView(i, len, indices);
    for(int j = 0; j < len; j++){
      indices[j] = A_->LCID(indices[j]);
    }
  }
}

int Ifpack_Euclid::Initialize(){
  comm_dh = GetMpiComm();
  MPI_Comm_size(comm_dh, &np_dh);
  MPI_Comm_rank(comm_dh, &myid_dh);
  Time_.ResetStartTime();
  if(mem_dh == NULL){
    Mem_dhCreate(&mem_dh);
  }
  if (tlog_dh == NULL) {
    TimeLog_dhCreate(&tlog_dh);
  }

  if (parser_dh == NULL) {
    Parser_dhCreate(&parser_dh);
  }
  Parser_dhInit(parser_dh, 0, NULL);
  //openLogfile_dh(0, NULL);
  Euclid_dhCreate(&eu);
  IsInitialized_=true;
  NumInitialize_ = NumInitialize_ + 1;
  InitializeTime_ = InitializeTime_ + Time_.ElapsedTime();
  return 0;
}

int Ifpack_Euclid::SetParameters(Teuchos::ParameterList& list){
  List_ = list;
  SetLevel_ = list.get("SetLevel", (int)1);
  SetBJ_ = list.get("SetBJ", (int)0);
  SetStats_ = list.get("SetStats", (int)0);
  SetMem_ = list.get("SetMem", (int)0);
  SetSparse_ = list.get("SetSparse", (double)0.0);
  SetRowScale_ = list.get("SetRowScale", (int)0);
  SetILUT_ = list.get("SetILUT", (double)0.0);
  return 0;
}

int Ifpack_Euclid::SetParameter(string name, int value){
  locale loc;
  for(size_t i = 0; i < name.length(); i++){
    name[i] = (char) tolower(name[i], loc);
  }
  if(name.compare("setlevel") == 0){
    SetLevel_ = value;
  } else if(name.compare("setbj") == 0){
    SetBJ_ = value;
  } else if(name.compare("setstats") == 0){
    SetStats_ = value;
  } else if(name.compare("setmem") == 0){
    SetMem_ = value;
  } else if(name.compare("setrowscale") == 0){
    SetRowScale_ = value;
  } else {
    cout << "\nThe string " << name << " is not an available option.\n";
    IFPACK_CHK_ERR(-1);
  }
  return 0;
}

int Ifpack_Euclid::SetParameter(string name, double value){
  locale loc;
  for(size_t i; i < name.length(); i++){
    name[i] = (char) tolower(name[i], loc);
  }
  if(name.compare("setsparse") == 0){
    SetSparse_ = value;
  } else if(name.compare("setilut") == 0){
    SetILUT_ = value;
  } else {
    cout << "\nThe string " << name << " is not an available option.\n";
    IFPACK_CHK_ERR(-1);
  }
  return 0;
}

int Ifpack_Euclid::Compute(){
  if(IsInitialized() == false){
    IFPACK_CHK_ERR(Initialize());
  }
  Time_.ResetStartTime();
  sprintf(Label_, "IFPACK_Euclid (level=%d, bj=%d, stats=%d, mem=%d, sparse = %f, rowscale = %d, ilut = %f)",
      SetLevel_, SetBJ_, SetStats_, SetMem_, SetSparse_, SetRowScale_, SetILUT_);

  // Set the parameters
  eu->level = SetLevel_;
  if(SetBJ_ != 0){
    strcpy("bj", eu->algo_par);
  }
  if(SetSparse_ != 0.0){
    eu->sparseTolA = SetSparse_;
  }
  if(SetRowScale_ != 0){
    eu->isScaled = true;
  }
  if(SetILUT_ != 0.0){
    eu->droptol = SetILUT_;
  }
  if(SetStats_ != 0 || SetMem_ != 0){
    eu->logging = true;
    Parser_dhInsert(parser_dh, "-eu_stats", "1");
  }
  eu->A = (void*) A_.get();
  eu->m = A_->NumMyRows();
  eu->n = A_->NumGlobalRows();
  Euclid_dhSetup(eu);
  IsComputed_ = true;
  NumCompute_ = NumCompute_ + 1;
  ComputeTime_ = ComputeTime_ + Time_.ElapsedTime();
  return 0;
}

const Epetra_Map& Ifpack_Euclid::OperatorDomainMap() const{
  return *MySimpleMap_;
}

const Epetra_Map& Ifpack_Euclid::OperatorRangeMap() const{
  return *MySimpleMap_;
}


int Ifpack_Euclid::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  if(IsComputed() == false){
    IFPACK_CHK_ERR(-1);
  }
  int NumVectors = X.NumVectors();
  if(NumVectors != Y.NumVectors()){
    IFPACK_CHK_ERR(-2);
  }
  Time_.ResetStartTime();
  for(int vecNum = 0; vecNum < NumVectors; vecNum++){
    CallEuclid(X[vecNum], Y[vecNum]);
  }
  Euclid_dhPrintTestData(eu, stdout);
  NumApplyInverse_ = NumApplyInverse_ + 1;
  ApplyInverseTime_ = ApplyInverseTime_ + Time_.ElapsedTime();
  return 0;
}

int Ifpack_Euclid::CallEuclid(double* x, double* y) const{
  Euclid_dhApply(eu, x, y);
  return 0;
}

ostream& operator << (ostream& os, const Ifpack_Euclid& A){
  if (!A.Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_Euclid: " << A.Label () << endl << endl;
    os << "Using " << A.Comm().NumProc() << " processors." << endl;
    os << "Global number of rows            = " << A.Matrix().NumGlobalRows() << endl;
    os << "Global number of nonzeros        = " << A.Matrix().NumGlobalNonzeros() << endl;
    os << "Condition number estimate = " << A.Condest() << endl;
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << A.NumInitialize()
       << "  " << std::setw(15) << A.InitializeTime()
       << "              0.0              0.0" << endl;
    os << "Compute()       "   << std::setw(5) << A.NumCompute()
       << "  " << std::setw(15) << A.ComputeTime()
       << "  " << std::setw(15) << 1.0e-6 * A.ComputeFlops();
    if (A.ComputeTime() != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * A.ComputeFlops() / A.ComputeTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << A.NumApplyInverse()
       << "  " << std::setw(15) << A.ApplyInverseTime()
       << "  " << std::setw(15) << 1.0e-6 * A.ApplyInverseFlops();
    if (A.ApplyInverseTime() != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * A.ApplyInverseFlops() / A.ApplyInverseTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }
  return os;
}


double Ifpack_Euclid::Condest(const Ifpack_CondestType CT, 
                             const int MaxIters,
                             const double Tol,
                             Epetra_RowMatrix* Matrix_in){
  if (!IsComputed()) // cannot compute right now
    return(-1.0);


  return(Condest_);
}

#endif // HAVE_EUCLID && HAVE_MPI
