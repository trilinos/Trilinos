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

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_IKLU.h"
#include "Ifpack_Condest.h"
#include "Ifpack_Utils.h"
#include "Ifpack_HashTable.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Util.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include <functional>
#include <algorithm>

using namespace Teuchos;

//==============================================================================
// FIXME: allocate Comm_ and Time_ the first Initialize() call
Ifpack_IKLU::Ifpack_IKLU(const Epetra_RowMatrix* A) :
  A_(*A),
  Comm_(A->Comm()),
  Condest_(-1.0),
  Relax_(0.),
  Athresh_(0.0),
  Rthresh_(1.0),
  LevelOfFill_(1.0),
  DropTolerance_(1e-12),
  IsInitialized_(false),
  IsComputed_(false),
  UseTranspose_(false),
  NumMyRows_(-1),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  Time_(Comm()),
  GlobalNonzeros_(0),
  csrA_(0),
  cssS_(0),
  csrnN_(0)
{
  // do nothing here..
}

//==============================================================================
Ifpack_IKLU::~Ifpack_IKLU()
{
  Destroy();
}

//==============================================================================
void Ifpack_IKLU::Destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
  if (csrA_)
    csr_spfree( csrA_ );
  if (cssS_)
    csr_sfree( cssS_ );
  if (csrnN_)
    csr_nfree( csrnN_ );
}

//==========================================================================
int Ifpack_IKLU::SetParameters(Teuchos::ParameterList& List)
{
  try 
  {
    LevelOfFill_ = List.get<double>("fact: ilut level-of-fill", LevelOfFill());
    if (LevelOfFill_ <= 0.0)
      IFPACK_CHK_ERR(-2); // must be greater than 0.0

    Athresh_ = List.get<double>("fact: absolute threshold", Athresh_);
    Rthresh_ = List.get<double>("fact: relative threshold", Rthresh_);
    Relax_ = List.get<double>("fact: relax value", Relax_);
    DropTolerance_ = List.get<double>("fact: drop tolerance", DropTolerance_);

    Label_ = "IFPACK IKLU (fill=" + Ifpack_toString(LevelOfFill())
      + ", relax=" + Ifpack_toString(RelaxValue())
      + ", athr=" + Ifpack_toString(AbsoluteThreshold())
      + ", rthr=" + Ifpack_toString(RelativeThreshold())
      + ", droptol=" + Ifpack_toString(DropTolerance())
      + ")";
    return(0);
  }
  catch (...)
  {
    cerr << "Caught an exception while parsing the parameter list" << endl;
    cerr << "This typically means that a parameter was set with the" << endl;
    cerr << "wrong type (for example, int instead of double). " << endl;
    cerr << "please check the documentation for the type required by each parameer." << endl;
    IFPACK_CHK_ERR(-1);
  }

  return(0);
}

//==========================================================================
int Ifpack_IKLU::Initialize()
{
  // delete previously allocated factorization
  Destroy();

  Time_.ResetStartTime();

  if (A_.Comm().NumProc() != 1) {
    cout << " There are too many processors !!! " << endl;
    cerr << "Ifpack_IKLU can be used with Comm().NumProc() == 1" << endl;
    cerr << "only. This class is a subdomain solver for Ifpack_AdditiveSchwarz," << endl;
    cerr << "and it is currently not meant to be used otherwise." << endl;
    exit(EXIT_FAILURE);
  }
  
  // check dimensions of input matrix only in serial
  if (Comm().NumProc() == 1 && Matrix().NumMyRows() != Matrix().NumMyCols())
    IFPACK_CHK_ERR(-2);
    
  NumMyRows_ = Matrix().NumMyRows();
  NumMyNonzeros_ = Matrix().NumMyNonzeros();

  int RowNnz, Length = Matrix().MaxNumEntries();
  std::vector<int>    RowIndices(Length);
  std::vector<double> RowValues(Length);

  //cout << "Processor " << Comm().MyPID() << " owns " << NumMyRows_ << " rows and has " << NumMyNonzeros_ << " nonzeros " << endl;
  // get general symbolic structure of the matrix
  csrA_ = csr_spalloc( NumMyRows_, NumMyRows_, NumMyNonzeros_, 1, 0 );

  // copy the symbolic structure into csrA_
  int count = 0;
  csrA_->p[0] = 0;
  for (int i = 0; i < NumMyRows_; ++i ) {

    IFPACK_CHK_ERR(A_.ExtractMyRowCopy(i,Length,RowNnz,
				       &RowValues[0],&RowIndices[0]));
    for (int j = 0 ; j < RowNnz ; ++j) {
      csrA_->j[count++] = RowIndices[j];
      //cout << "Row = " << i << ", Column = " << RowIndices[j] << ", Value = " << RowValues[j] << endl;
    }
    csrA_->p[i+1] = csrA_->p[i] + RowNnz;
  }

  // Perform symbolic analysis on the current matrix structure
  int order = 1;
  cssS_ = csr_sqr( order, csrA_ );
  for (int i = 0; i < NumMyRows_; ++i ) {
     cout << "AMD Perm (from inside KLU) [" << i << "] = " << cssS_->q[i] << endl;
  }

  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_.ElapsedTime();

  return(0);
}

//==========================================================================
class Ifpack_AbsComp 
{
 public:
  inline bool operator()(const double& x, const double& y) 
  {
    return(IFPACK_ABS(x) > IFPACK_ABS(y));
  }
};

//==========================================================================

int Ifpack_IKLU::Compute() 
{
  if (!IsInitialized()) 
    IFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  NumMyRows_ = A_.NumMyRows();
  int Length = A_.MaxNumEntries();

  bool distributed = (Comm().NumProc() > 1)?true:false;
#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  if (distributed)
  {
    SerialComm_ = rcp(new Epetra_SerialComm);
    SerialMap_ = rcp(new Epetra_Map(NumMyRows_, 0, *SerialComm_));
    assert (SerialComm_.get() != 0);
    assert (SerialMap_.get() != 0);
  }
  else
    SerialMap_ = rcp(const_cast<Epetra_Map*>(&A_.RowMatrixRowMap()), false);
#endif

  int RowNnz;
  std::vector<int>    RowIndices(Length);
  std::vector<double> RowValues(Length);

  // copy the values from A_ into csrA_
  int count = 0;
  for (int i = 0; i < NumMyRows_; ++i ) {

    IFPACK_CHK_ERR(A_.ExtractMyRowCopy(i,Length,RowNnz,
				       &RowValues[0],&RowIndices[0]));
    // make sure each row has the same number of nonzeros
    if (RowNnz != (csrA_->p[i+1]-csrA_->p[i])) {
      cout << "The number of nonzeros for this row does not math the expected number of nonzeros!!!" << endl;
    }
    for (int j = 0 ; j < RowNnz ; ++j) {
      
      csrA_->x[count++] = RowValues[j];
      //cout << "Row = " << i << ", Column = " << RowIndices[j] << ", Value = " << RowValues[j] << endl;
    }
  }
  
  // compute the lu factors
  double tol = 0.1;
  csrnN_ = csr_lu( &*csrA_, &*cssS_, tol );

  // Create L and U as a view of the information stored in csrnN_->L and csrnN_->U
  csr* L_tmp = csrnN_->L;
  csr* U_tmp = csrnN_->U;
  std::vector<int> numEntriesL( NumMyRows_ ), numEntriesU( NumMyRows_ );
  for (int i=0; i < NumMyRows_; ++i) {
    numEntriesL[i] = ( L_tmp->p[i+1] - L_tmp->p[i] );
    numEntriesU[i] = ( U_tmp->p[i+1] - U_tmp->p[i] );
  }
  L_ = rcp(new Epetra_CrsMatrix(View, *SerialMap_, &numEntriesL[0]));
  U_ = rcp(new Epetra_CrsMatrix(View, *SerialMap_, &numEntriesU[0]));

  // Insert the values into L and U
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(SerialMap_->GlobalIndicesInt()) {
      for (int i=0; i < NumMyRows_; ++i) {
        L_->InsertGlobalValues( i, numEntriesL[i], &(L_tmp->x[L_tmp->p[i]]), &(L_tmp->j[L_tmp->p[i]]) );
        U_->InsertGlobalValues( i, numEntriesU[i], &(U_tmp->x[U_tmp->p[i]]), &(U_tmp->j[U_tmp->p[i]]) );
      }
    }
    else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    if(SerialMap_->GlobalIndicesLongLong()) {

      const int MaxNumEntries_L_U = std::max(L_->MaxNumEntries(), U_->MaxNumEntries());
	  std::vector<long long> entries(MaxNumEntries_L_U);

      for (int i=0; i < NumMyRows_; ++i) {
        std::copy(&(L_tmp->j[L_tmp->p[i]]), &(L_tmp->j[L_tmp->p[i]]) + numEntriesL[i], entries.begin());
        L_->InsertGlobalValues( i, numEntriesL[i], &(L_tmp->x[L_tmp->p[i]]), &(entries[0]) );

        std::copy(&(U_tmp->j[U_tmp->p[i]]), &(U_tmp->j[U_tmp->p[i]]) + numEntriesU[i], entries.begin());
        U_->InsertGlobalValues( i, numEntriesU[i], &(U_tmp->x[U_tmp->p[i]]), &(entries[0]) );
      }
    }
    else
#endif
      throw "Ifpack_IKLU::Compute: GlobalIndices type unknown for SerialMap_";

  IFPACK_CHK_ERR(L_->FillComplete());
  IFPACK_CHK_ERR(U_->FillComplete());

  long long MyNonzeros = L_->NumGlobalNonzeros64() + U_->NumGlobalNonzeros64();
  Comm().SumAll(&MyNonzeros, &GlobalNonzeros_, 1);

  IsComputed_ = true;

  ++NumCompute_;
  ComputeTime_ += Time_.ElapsedTime();

  return(0);

}
  
//=============================================================================
int Ifpack_IKLU::ApplyInverse(const Epetra_MultiVector& X, 
			     Epetra_MultiVector& Y) const
{
  if (!IsComputed())
    IFPACK_CHK_ERR(-2); // compute preconditioner first

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK_CHK_ERR(-3); // Return error: X and Y not the same size

  Time_.ResetStartTime();

  // NOTE: L_ and U_ are based on SerialMap_, while Xcopy is based
  // on A.Map()... which are in general different. However, Solve()
  // does not seem to care... which is fine with me.
  //
  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy and apply permutation.
  std::vector<int> invq( NumMyRows_ );

  for (int i=0; i<NumMyRows_; ++i ) {
    csrnN_->perm[ csrnN_->pinv[i] ] = i;
    invq[ cssS_->q[i] ] = i;
  }

  Teuchos::RefCountPtr<Epetra_MultiVector> Xcopy = Teuchos::rcp( new Epetra_MultiVector(X.Map(),X.NumVectors()), false );
  Teuchos::RefCountPtr<Epetra_MultiVector> Ytemp = Teuchos::rcp( new Epetra_MultiVector(Y.Map(),Y.NumVectors()) );

  for (int i=0; i<NumMyRows_; ++i) {
    for (int j=0; j<X.NumVectors(); ++j) {
      Xcopy->ReplaceMyValue( invq[i], j, (*X(j))[i] );
    }
  }

  if (!UseTranspose_)
  {
    // solves LU Y = X 
    IFPACK_CHK_ERR(L_->Solve(false,false,false,*Xcopy,*Ytemp));
    IFPACK_CHK_ERR(U_->Solve(true,false,false,*Ytemp,*Ytemp));
  }
  else
  {
    // solves U(trans) L(trans) Y = X
    IFPACK_CHK_ERR(U_->Solve(true,true,false,*Xcopy,*Ytemp));
    IFPACK_CHK_ERR(L_->Solve(false,true,false,*Ytemp,*Ytemp));
  }

  // Reverse the permutation.
  for (int i=0; i<NumMyRows_; ++i) {
    for (int j=0; j<Y.NumVectors(); ++j) {
      Y.ReplaceMyValue( csrnN_->perm[i], j, (*(*Ytemp)(j))[i] );
    }
  }

  ++NumApplyInverse_;
#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += X.NumVectors() * 2 * GlobalNonzeros_;
#endif
  ApplyInverseTime_ += Time_.ElapsedTime();

  return(0);

}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_IKLU::Apply(const Epetra_MultiVector& X, 
		      Epetra_MultiVector& Y) const 
{

  return(-98);
}

//=============================================================================
double Ifpack_IKLU::Condest(const Ifpack_CondestType CT, 
                            const int MaxIters, const double Tol,
			    Epetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // NOTE: this is computing the *local* condest
  if (Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//=============================================================================
std::ostream&
Ifpack_IKLU::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_IKLU: " << Label() << endl << endl;
    os << "Level-of-fill      = " << LevelOfFill() << endl;
    os << "Absolute threshold = " << AbsoluteThreshold() << endl;
    os << "Relative threshold = " << RelativeThreshold() << endl;
    os << "Relax value        = " << RelaxValue() << endl;
    os << "Condition number estimate       = " << Condest() << endl;
    os << "Global number of rows           = " << A_.NumGlobalRows64() << endl;
    if (IsComputed_) {
      os << "Number of nonzeros in A         = " << A_.NumGlobalNonzeros64() << endl;
      os << "Number of nonzeros in L + U     = " << NumGlobalNonzeros64() 
         << " ( = " << 100.0 * NumGlobalNonzeros64() / A_.NumGlobalNonzeros64() 
         << " % of A)" << endl;
      os << "nonzeros / rows                 = " 
        << 1.0 * NumGlobalNonzeros64() / U_->NumGlobalRows64() << endl;
    }
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize() 
       << "  " << std::setw(15) << InitializeTime() 
       << "               0.0            0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute() 
       << "  " << std::setw(15) << ComputeTime()
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops(); 
    if (ComputeTime() != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops() / ComputeTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse() 
       << "  " << std::setw(15) << ApplyInverseTime()
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops();
    if (ApplyInverseTime() != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops() / ApplyInverseTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  return(os);
}
