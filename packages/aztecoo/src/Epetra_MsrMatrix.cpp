
/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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

#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_MsrMatrix.h"
//==============================================================================
Epetra_MsrMatrix::Epetra_MsrMatrix(int * proc_config, AZ_MATRIX * a_mat)
  : Epetra_Object("Epetra::MsrMatrix"),
    Amat_(a_mat),
    proc_config_(proc_config),
    Values_(0),
    Indices_(0),
    MaxNumEntries_(-1),
    ImportVector_(0),
    NormInf_(-1.0),
    NormOne_(-1.0)
{
#ifdef AZTEC_MPI
  MPI_Comm * mpicomm = (MPI_Comm * ) AZ_get_comm(proc_config);
  Comm_ = new Epetra_MpiComm(*mpicomm);
#else
  Comm_ = new Epetra_SerialComm();
#endif  
  if (a_mat->data_org[AZ_matrix_type]!=AZ_MSR_MATRIX)
    throw Comm_->ReportError("AZ_matrix_type must be AZ_MSR_MATRIX", -1);
  int * bindx = a_mat->bindx;
  NumMyRows_ = a_mat->data_org[AZ_N_internal] + a_mat->data_org[AZ_N_border];
  int NumExternal = a_mat->data_org[AZ_N_external];
  NumMyCols_ = NumMyRows_ + NumExternal;
  NumMyNonzeros_ = bindx[NumMyRows_] - bindx[0] + NumMyRows_;
  //Comm_->SumAll(&NumMyNonzeros_, &NumGlobalNonzeros_, 1);
  long long NumMyNonzerosLL_ = (long long) NumMyNonzeros_;
  Comm_->SumAll(&NumMyNonzerosLL_, &NumGlobalNonzeros_, 1);

  int * MyGlobalElements = a_mat->update;
  if (MyGlobalElements==0) 
    throw Comm_->ReportError("Aztec matrix has no update list: Check if AZ_Transform was called.", -2);

  DomainMap_ = new Epetra_Map(-1, NumMyRows_, MyGlobalElements, 0, *Comm_);

  double * dbleColGIDs = new double[NumMyCols_];
  int * ColGIDs = new int[NumMyCols_];
  for (int i=0; i<NumMyRows_; i++) dbleColGIDs[i] = (double) MyGlobalElements[i];
  AZ_exchange_bdry(dbleColGIDs, a_mat->data_org, proc_config);
  {for (int i=0; i<NumMyCols_; i++) ColGIDs[i] = (int) dbleColGIDs[i];}

  ColMap_ = new Epetra_Map(-1, NumMyCols_, ColGIDs, 0, *Comm_);

  Importer_ = new Epetra_Import(*ColMap_, *DomainMap_);
  
  delete [] dbleColGIDs;
  delete [] ColGIDs;
}
//==============================================================================
Epetra_MsrMatrix::~Epetra_MsrMatrix(){
  if (ImportVector_!=0) delete ImportVector_;
  delete Importer_;
  delete ColMap_;
  delete DomainMap_;
  delete Comm_;

  if (Values_!=0) {delete [] Values_; Values_=0;}
  if (Indices_!=0) {delete [] Indices_; Indices_=0;}
}
//==========================================================================
int Epetra_MsrMatrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, double * Values,
					 int * Indices) const 
{
  int tmpNumEntries;
  int err = AZ_MSR_getrow(Indices, Values, &tmpNumEntries, Amat_, 1, &Row, Length);
  NumEntries = tmpNumEntries;
  return(err-1);
}
//==========================================================================
int Epetra_MsrMatrix::NumMyRowEntries(int Row, int & NumEntries) const 
{

  if (RowMatrixRowMap().MyLID(Row)) NumEntries = Amat_->bindx[Row+1] - Amat_->bindx[Row] + 1;
  else {
    EPETRA_CHK_ERR(-1);
  }
  return(0);
}
//==============================================================================
int Epetra_MsrMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const {
	
  int iend = NumMyDiagonals();
  for (int i=0; i<iend; i++) Diagonal[i] = Amat_->val[i];
  return(0);
}
//=============================================================================
int Epetra_MsrMatrix::Multiply(bool TransA,
                               const Epetra_MultiVector& X,
                               Epetra_MultiVector& Y) const
{
  (void)TransA;
  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) EPETRA_CHK_ERR(-1);  // X and Y must have same number of vectors

  double ** xptrs;
  double ** yptrs;
  X.ExtractView(&xptrs);
  Y.ExtractView(&yptrs);
  if (RowMatrixImporter()!=0) {
    if (ImportVector_!=0) {
      if (ImportVector_->NumVectors()!=NumVectors) { delete ImportVector_; ImportVector_= 0;}
    }
    if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(RowMatrixColMap(),NumVectors); // Create import vector if needed
    ImportVector_->Import(X, *RowMatrixImporter(), Insert);
    ImportVector_->ExtractView(&xptrs);
  }
  for (int i=0; i<NumVectors; i++)
    Amat_->matvec(xptrs[i], yptrs[i], Amat_, proc_config_);
  
  double flops = NumGlobalNonzeros();
  flops *= 2.0;
  flops *= (double) NumVectors;
  UpdateFlops(flops);
  return(0);
}

//=============================================================================
int Epetra_MsrMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, 
                            const Epetra_MultiVector& X,
                            Epetra_MultiVector& Y) const
{
  (void)Upper;
  (void)Trans;
  (void)UnitDiagonal;
  (void)X;
  (void)Y;
  return(-1); // Not implemented
}

//=============================================================================
// Utility routine to get the specified row and put it into Values_ and Indices_ arrays
int Epetra_MsrMatrix::MaxNumEntries() const {
  int NumEntries;

  if (MaxNumEntries_==-1) {
    for (int i=0; i < NumMyRows_; i++) {
      NumMyRowEntries(i, NumEntries);
      if (NumEntries>MaxNumEntries_) MaxNumEntries_ = NumEntries;
    }
  }
  return(MaxNumEntries_);
}

//=============================================================================
// Utility routine to get the specified row and put it into Values_ and Indices_ arrays
int Epetra_MsrMatrix::GetRow(int Row) const {

  int NumEntries;
  int maxNumEntries = Epetra_MsrMatrix::MaxNumEntries();

  if (maxNumEntries>0) {
    if (Values_==0) Values_ = new double[maxNumEntries];
    if (Indices_==0) Indices_ = new int[maxNumEntries];
  }
  Epetra_MsrMatrix::ExtractMyRowCopy(Row, maxNumEntries, NumEntries, Values_, Indices_);
  
  return(NumEntries);
}

//=============================================================================
int Epetra_MsrMatrix::InvRowSums(Epetra_Vector& x) const {
//
// Put inverse of the sum of absolute values of the ith row of A in x[i].
//

  if (!OperatorRangeMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the range of A

  int ierr = 0;
  int i, j;
  for (i=0; i < NumMyRows_; i++) {
    int NumEntries = GetRow(i); // Copies ith row of matrix into Values_ and Indices_
    double scale = 0.0;
    for (j=0; j < NumEntries; j++) scale += fabs(Values_[j]);
    if (scale<Epetra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      x[i] = Epetra_MaxDouble;
    }
    else
      x[i] = 1.0/scale;
  }
  UpdateFlops(NumGlobalNonzeros());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
//=============================================================================
int Epetra_MsrMatrix::InvColSums(Epetra_Vector& x) const {
//
// Put inverse of the sum of absolute values of the jth column of A in x[j].
//

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  if (!OperatorDomainMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A
  

  Epetra_Vector * xp = 0;
  Epetra_Vector * x_tmp = 0;
  

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if (RowMatrixImporter()!=0) {
    x_tmp = new Epetra_Vector(RowMatrixColMap()); // Create import vector if needed
    xp = x_tmp;
  }
  int ierr = 0;
  int i, j;

  for (i=0; i < NumMyCols_; i++) (*xp)[i] = 0.0;

  for (i=0; i < NumMyRows_; i++) {
    int NumEntries = GetRow(i);// Copies ith row of matrix into Values_ and Indices_
    for (j=0; j < NumEntries; j++) (*xp)[Indices_[j]] += fabs(Values_[j]);
  }

  if (RowMatrixImporter()!=0){
    x.Export(*x_tmp, *RowMatrixImporter(), Add); // Fill x with Values from import vector
    delete x_tmp;
    xp = &x;
  }
  // Invert values, don't allow them to get too large
  for (i=0; i < NumMyRows_; i++) {
    double scale = (*xp)[i];
    if (scale<Epetra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      (*xp)[i] = Epetra_MaxDouble;
    }
    else
      (*xp)[i] = 1.0/scale;
  }
  UpdateFlops(NumGlobalNonzeros());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_MsrMatrix::LeftScale(const Epetra_Vector& x) {
//
// This function scales the ith row of A by x[i].
//
  // For this method, we have no choice but to work with the UGLY MSR data structures.

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  if (!OperatorRangeMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the range of A

  int i, j;
  int * bindx = Amat_->bindx;
  double * val = Amat_->val;


  for (i=0; i < NumMyRows_; i++) {
    
    int NumEntries = bindx[i+1] - bindx[i];
    double scale = x[i];
    val[i] *= scale;
    double * Values = val + bindx[i];
    for (j=0; j < NumEntries; j++)  Values[j] *= scale;
  }
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Epetra_MsrMatrix::RightScale(const Epetra_Vector& x) {
//
// This function scales the jth row of A by x[j].
//

  // For this method, we have no choice but to work with the UGLY MSR data structures.

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  if (!OperatorDomainMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A

  int * bindx = Amat_->bindx;
  double * val = Amat_->val;
  Epetra_Vector * xp = 0;
  Epetra_Vector * x_tmp = 0;

  // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
  if (RowMatrixImporter()!=0) {
    x_tmp = new Epetra_Vector(RowMatrixColMap()); // Create import vector if needed
    x_tmp->Import(x,*RowMatrixImporter(), Insert); // x_tmp will have all the values we need
    xp = x_tmp;
  }

  int i, j;

  for (i=0; i < NumMyRows_; i++) {
    int NumEntries = bindx[i+1] - bindx[i];
    double scale = (*xp)[i];
    val[i] *= scale;
    double * Values = val + bindx[i];
    int * Indices = bindx + bindx[i];
    for (j=0; j < NumEntries; j++)  Values[j] *=  (*xp)[Indices[j]];
  }
  if (x_tmp!=0) delete x_tmp;
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}

//=============================================================================
double Epetra_MsrMatrix::NormInf() const {

  if (NormInf_>-1.0) return(NormInf_);

  double Local_NormInf = 0.0;
  for (int i=0; i < NumMyRows_; i++) {
    int NumEntries = GetRow(i);
    double sum = 0.0;
    for (int j=0; j < NumEntries; j++) sum += fabs(Values_[j]);
    
    Local_NormInf = EPETRA_MAX(Local_NormInf, sum);
  }
  Comm().MaxAll(&Local_NormInf, &NormInf_, 1);
  UpdateFlops(NumGlobalNonzeros());
  return(NormInf_);
}
//=============================================================================
double Epetra_MsrMatrix::NormOne() const {

  if (NormOne_>-1.0) return(NormOne_);

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  Epetra_Vector * x = new Epetra_Vector(RowMatrixRowMap()); // Need temp vector for column sums
  
  Epetra_Vector * xp = 0;
  Epetra_Vector * x_tmp = 0;
  

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if (RowMatrixImporter()!=0) {
    x_tmp = new Epetra_Vector(RowMatrixColMap()); // Create temporary import vector if needed
    xp = x_tmp;
  }
  int i, j;

  for (i=0; i < NumMyCols_; i++) (*xp)[i] = 0.0;

  for (i=0; i < NumMyRows_; i++) {
    int NumEntries = GetRow(i);
    for (j=0; j < NumEntries; j++) (*xp)[Indices_[j]] += fabs(Values_[j]);
  }
  if (RowMatrixImporter()!=0) x->Export(*x_tmp, *RowMatrixImporter(), Add); // Fill x with Values from temp vector
  x->MaxValue(&NormOne_); // Find max
  if (x_tmp!=0) delete x_tmp;
  delete x;
  UpdateFlops(NumGlobalNonzeros());
  return(NormOne_);
}
//=============================================================================
void Epetra_MsrMatrix::Print(ostream& os) const {
  int MyPID = RowMatrixRowMap().Comm().MyPID();
  int NumProc = RowMatrixRowMap().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
/*      const Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
      const Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
      const int             oldp = os.precision(12); */
      if (MyPID==0) {
	os <<  "\nNumber of Global Rows        = "; os << NumGlobalRows(); os << endl;
	os <<    "Number of Global Cols        = "; os << NumGlobalCols(); os << endl;
	os <<    "Number of Global Diagonals   = "; os << NumGlobalDiagonals(); os << endl;
	os <<    "Number of Global Nonzeros    = "; os << NumGlobalNonzeros(); os << endl;
	if (LowerTriangular()) os <<    " ** Matrix is Lower Triangular **"; os << endl;
	if (UpperTriangular()) os <<    " ** Matrix is Upper Triangular **"; os << endl;
      }

      os <<  "\nNumber of My Rows        = "; os << NumMyRows(); os << endl;
      os <<    "Number of My Cols        = "; os << NumMyCols(); os << endl;
      os <<    "Number of My Diagonals   = "; os << NumMyDiagonals(); os << endl;
      os <<    "Number of My Nonzeros    = "; os << NumMyNonzeros(); os << endl; os << endl;

      os << flush;
      
      // Reset os flags
      
/*      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp); */
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }

  {for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int i, j;

      if (MyPID==0) {
	os.width(8);
	os <<  "   Processor ";
	os.width(10);
	os <<  "   Row Index ";
	os.width(10);
	os <<  "   Col Index ";
	os.width(20);
	os <<  "   Value     ";
	os << endl;
      }
      for (i=0; i<NumMyRows_; i++) {
	int Row = RowMatrixRowMap().GID(i); // Get global row number
	int NumEntries = GetRow(i); // ith row is now in Values_ and Indices_
	
	for (j = 0; j < NumEntries ; j++) {   
	  os.width(8);
	  os <<  MyPID ; os << "    ";	
	  os.width(10);
	  os <<  Row ; os << "    ";	
	  os.width(10);
	  os <<  RowMatrixColMap().GID(Indices_[j]); os << "    ";
	  os.width(20);
	  os <<  Values_[j]; os << "    ";
	  os << endl;
	}
      }

      
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    RowMatrixRowMap().Comm().Barrier();
    RowMatrixRowMap().Comm().Barrier();
    RowMatrixRowMap().Comm().Barrier();
  }}

  return;
}
