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
#include "Ifpack_ILUT.h"
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
Ifpack_ILUT::Ifpack_ILUT(const Epetra_RowMatrix* A) :
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
  GlobalNonzeros_(0)
{
  // do nothing here..
}

//==============================================================================
Ifpack_ILUT::~Ifpack_ILUT()
{
  Destroy();
}

//==============================================================================
void Ifpack_ILUT::Destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
}

//==========================================================================
int Ifpack_ILUT::SetParameters(Teuchos::ParameterList& List)
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

    Label_ = "IFPACK ILUT (fill=" + Ifpack_toString(LevelOfFill())
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
int Ifpack_ILUT::Initialize()
{
  // delete previously allocated factorization
  Destroy();

  Time_.ResetStartTime();

  // check only in serial
  if (Comm().NumProc() == 1 && Matrix().NumMyRows() != Matrix().NumMyCols())
    IFPACK_CHK_ERR(-2);
    
  NumMyRows_ = Matrix().NumMyRows();

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
// MS // completely rewritten the algorithm on 25-Jan-05, using more STL
// MS // and in a nicer and cleaner way. Also, it is more efficient.
// MS // WARNING: Still not fully tested!
template<typename int_type>
int Ifpack_ILUT::TCompute() 
{
  int Length = A_.MaxNumEntries();
  std::vector<double> RowValuesL(Length);
  std::vector<int>    RowIndicesU(Length);
  std::vector<int_type> RowIndicesU_LL;
  RowIndicesU_LL.resize(Length);
  std::vector<double> RowValuesU(Length);
  
  int RowNnzU;

  L_ = rcp(new Epetra_CrsMatrix(Copy, *SerialMap_, 0));
  U_ = rcp(new Epetra_CrsMatrix(Copy, *SerialMap_, 0));

  if ((L_.get() == 0) || (U_.get() == 0))
    IFPACK_CHK_ERR(-5); // memory allocation error

  // insert first row in U_ and L_
  IFPACK_CHK_ERR(A_.ExtractMyRowCopy(0,Length,RowNnzU,
                                     &RowValuesU[0],&RowIndicesU[0]));

  bool distributed = (Comm().NumProc() > 1)?true:false;

  if (distributed)
  {
    int count = 0;
    for (int i = 0 ;i < RowNnzU ; ++i) 
    {
      if (RowIndicesU[i] < NumMyRows_){
        RowIndicesU[count] = RowIndicesU[i];
        RowValuesU[count] = RowValuesU[i];
        ++count;
      }
      else
        continue;
    }
    RowNnzU = count;
  }

  // modify diagonal
  for (int i = 0 ;i < RowNnzU ; ++i) {
    if (RowIndicesU[i] == 0) {
      double& v = RowValuesU[i];
      v = AbsoluteThreshold() * EPETRA_SGN(v) + RelativeThreshold() * v;
      break;
    }
  }
  
  std::copy(&(RowIndicesU[0]), &(RowIndicesU[0]) + RowNnzU, RowIndicesU_LL.begin());
  IFPACK_CHK_ERR(U_->InsertGlobalValues(0,RowNnzU,&(RowValuesU[0]),
                                        &(RowIndicesU_LL[0])));
   // FIXME: DOES IT WORK IN PARALLEL ??
  RowValuesU[0] = 1.0;
  RowIndicesU[0] = 0;
  int_type RowIndicesU_0_templ = RowIndicesU[0];
  IFPACK_CHK_ERR(L_->InsertGlobalValues(0,1,&(RowValuesU[0]),
                                        &RowIndicesU_0_templ));

  int max_keys =  (int) 10 * A_.MaxNumEntries() * LevelOfFill() ;
  Ifpack_HashTable SingleRowU(max_keys , 1);
  Ifpack_HashTable SingleRowL(max_keys , 1);

  int hash_size = SingleRowU.getRecommendedHashSize(max_keys) ;
  std::vector<int> keys;      keys.reserve(hash_size * 10);
  std::vector<double> values; values.reserve(hash_size * 10);
  std::vector<double> AbsRow; AbsRow.reserve(hash_size * 10);

  // =================== //
  // start factorization //
  // =================== //
  
#ifdef IFPACK_FLOPCOUNTERS
  double this_proc_flops = 0.0;
#endif

  for (int row_i = 1 ; row_i < NumMyRows_ ; ++row_i) 
  {
    // get row `row_i' of the matrix, store in U pointers
    IFPACK_CHK_ERR(A_.ExtractMyRowCopy(row_i,Length,RowNnzU,
                                       &RowValuesU[0],&RowIndicesU[0]));

    if (distributed)
    {
      int count = 0;
      for (int i = 0 ;i < RowNnzU ; ++i) 
      {
        if (RowIndicesU[i] < NumMyRows_){
          RowIndicesU[count] = RowIndicesU[i];
          RowValuesU[count] = RowValuesU[i];
          ++count;
        }
        else
          continue;
      }
      RowNnzU = count;
    }

    int NnzLower = 0;
    int NnzUpper = 0;

    for (int i = 0 ;i < RowNnzU ; ++i) {
      if (RowIndicesU[i] < row_i)
        NnzLower++;
      else if (RowIndicesU[i] == row_i) {
        // add threshold
        NnzUpper++;
        double& v = RowValuesU[i];
        v = AbsoluteThreshold() * EPETRA_SGN(v) + RelativeThreshold() * v;
      }
      else
        NnzUpper++;
    }

    int FillL = (int)(LevelOfFill() * NnzLower);
    int FillU = (int)(LevelOfFill() * NnzUpper);
    if (FillL == 0) FillL = 1;
    if (FillU == 0) FillU = 1;

    // convert line `row_i' into STL map for fast access
    SingleRowU.reset();

    for (int i = 0 ; i < RowNnzU ; ++i) {
        SingleRowU.set(RowIndicesU[i], RowValuesU[i]);
    }

    // for the multipliers
    SingleRowL.reset();

    int start_col = NumMyRows_;
    for (int i = 0 ; i < RowNnzU ; ++i)
      start_col = EPETRA_MIN(start_col, RowIndicesU[i]);

#ifdef IFPACK_FLOPCOUNTERS
    int flops = 0;
#endif
  
    for (int col_k = start_col ; col_k < row_i ; ++col_k) {

      int_type* ColIndicesK;
      double* ColValuesK;
      int     ColNnzK;

      double xxx = SingleRowU.get(col_k);
      // This factorization is too "relaxed". Leaving it as it is, as Tifpack
      // does it correctly.
      if (IFPACK_ABS(xxx) > DropTolerance()) {
          IFPACK_CHK_ERR(U_->ExtractGlobalRowView(col_k, ColNnzK, ColValuesK, 
                                                  ColIndicesK));

          // FIXME: can keep trace of diagonals
          double DiagonalValueK = 0.0;
          for (int i = 0 ; i < ColNnzK ; ++i) {
            if (ColIndicesK[i] == col_k) {
              DiagonalValueK = ColValuesK[i];
              break;
            }
          }

          // FIXME: this should never happen!
          if (DiagonalValueK == 0.0)
            DiagonalValueK = AbsoluteThreshold();

          double yyy = xxx / DiagonalValueK ;
          SingleRowL.set(col_k, yyy);
#ifdef IFPACK_FLOPCOUNTERS
          ++flops;
#endif

          for (int j = 0 ; yyy != 0.0 && j < ColNnzK ; ++j)
          {
            int_type col_j = ColIndicesK[j];

            if (col_j < col_k) continue;

            SingleRowU.set(col_j, -yyy * ColValuesK[j], true);
#ifdef IFPACK_FLOPCOUNTERS
            flops += 2;
#endif
          }
      }
    }

#ifdef IFPACK_FLOPCOUNTERS
    this_proc_flops += 1.0 * flops;
#endif

    double cutoff = DropTolerance();
    double DiscardedElements = 0.0;
    int count;

    // drop elements to satisfy LevelOfFill(), start with L
    count = 0;
    int sizeL = SingleRowL.getNumEntries();
    keys.resize(sizeL);
    values.resize(sizeL);

    AbsRow.resize(sizeL);

    SingleRowL.arrayify(
      keys.size() ? &keys[0] : 0,
      values.size() ? &values[0] : 0
      );
    for (int i = 0; i < sizeL; ++i)
      if (IFPACK_ABS(values[i]) > DropTolerance()) {
        AbsRow[count++] = IFPACK_ABS(values[i]);
      }

    if (count > FillL) {
      nth_element(AbsRow.begin(), AbsRow.begin() + FillL, AbsRow.begin() + count, 
                  std::greater<double>());
      cutoff = AbsRow[FillL];
    }

    for (int i = 0; i < sizeL; ++i) {
      if (IFPACK_ABS(values[i]) >= cutoff) {
        int_type key_templ = keys[i];
        IFPACK_CHK_ERR(L_->InsertGlobalValues(row_i,1, &values[i], &key_templ));
      }
      else
        DiscardedElements += values[i];
    }

    // FIXME: DOES IT WORK IN PARALLEL ???
    // add 1 to the diagonal
    double dtmp = 1.0;
	int_type row_i_templ = row_i;
    IFPACK_CHK_ERR(L_->InsertGlobalValues(row_i,1, &dtmp, &row_i_templ));

    // same business with U_
    count = 0;
    int sizeU = SingleRowU.getNumEntries();
    AbsRow.resize(sizeU + 1);
    keys.resize(sizeU + 1);
    values.resize(sizeU + 1);

    SingleRowU.arrayify(&keys[0], &values[0]);

    for (int i = 0; i < sizeU; ++i)
      if (keys[i] >= row_i && IFPACK_ABS(values[i]) > DropTolerance())
      {
        AbsRow[count++] = IFPACK_ABS(values[i]);
      }

    if (count > FillU) {
      nth_element(AbsRow.begin(), AbsRow.begin() + FillU, AbsRow.begin() + count, 
                  std::greater<double>());
      cutoff = AbsRow[FillU];
    }

    // sets the factors in U_
    for (int i = 0; i < sizeU; ++i) 
    {
      int col = keys[i];
      double val = values[i];

      if (col >= row_i) {
        if (IFPACK_ABS(val) >= cutoff || row_i == col) {
          int_type col_templ = col;
          IFPACK_CHK_ERR(U_->InsertGlobalValues(row_i,1, &val, &col_templ));
        }
        else
          DiscardedElements += val;
      }
    }

    // FIXME: not so sure of that!
    if (RelaxValue() != 0.0) {
      DiscardedElements *= RelaxValue();
      IFPACK_CHK_ERR(U_->InsertGlobalValues(row_i,1, &DiscardedElements,
                                            &row_i_templ));
    }
  }

#ifdef IFPACK_FLOPCOUNTERS
  double tf;
  Comm().SumAll(&this_proc_flops, &tf, 1);
  ComputeFlops_ += tf;
#endif

  IFPACK_CHK_ERR(L_->FillComplete());
  IFPACK_CHK_ERR(U_->FillComplete());

#if 0
  // to check the complete factorization
  Epetra_Vector LHS(A_.RowMatrixRowMap());
  Epetra_Vector RHS1(A_.RowMatrixRowMap());
  Epetra_Vector RHS2(A_.RowMatrixRowMap());
  Epetra_Vector RHS3(A_.RowMatrixRowMap());
  LHS.Random();

  cout << "A = " << A_.NumGlobalNonzeros() << endl;
  cout << "L = " << L_->NumGlobalNonzeros() << endl;
  cout << "U = " << U_->NumGlobalNonzeros() << endl;

  A_.Multiply(false,LHS,RHS1);
  U_->Multiply(false,LHS,RHS2);
  L_->Multiply(false,RHS2,RHS3);

  RHS1.Update(-1.0, RHS3, 1.0);
  double Norm;
  RHS1.Norm2(&Norm);
#endif

  long long MyNonzeros = L_->NumGlobalNonzeros64() + U_->NumGlobalNonzeros64();
  Comm().SumAll(&MyNonzeros, &GlobalNonzeros_, 1);

  IsComputed_ = true;

  ++NumCompute_;
  ComputeTime_ += Time_.ElapsedTime();

  return(0);

}

int Ifpack_ILUT::Compute() {
  if (!IsInitialized()) 
    IFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  NumMyRows_ = A_.NumMyRows();
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

  int ret_val = 1;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(SerialMap_->GlobalIndicesInt()) {
    ret_val = TCompute<int>();
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(SerialMap_->GlobalIndicesLongLong()) {
    ret_val = TCompute<long long>();
  }
  else
#endif
    throw "Ifpack_ILUT::Compute: GlobalIndices type unknown";

  return ret_val;
}

//=============================================================================
int Ifpack_ILUT::ApplyInverse(const Epetra_MultiVector& X, 
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
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr<const Epetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  if (!UseTranspose_)
  {
    // solves LU Y = X 
    IFPACK_CHK_ERR(L_->Solve(false,false,false,*Xcopy,Y));
    IFPACK_CHK_ERR(U_->Solve(true,false,false,Y,Y));
  }
  else
  {
    // solves U(trans) L(trans) Y = X
    IFPACK_CHK_ERR(U_->Solve(true,true,false,*Xcopy,Y));
    IFPACK_CHK_ERR(L_->Solve(false,true,false,Y,Y));
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
int Ifpack_ILUT::Apply(const Epetra_MultiVector& X, 
		      Epetra_MultiVector& Y) const 
{
  return(-98);
}

//=============================================================================
double Ifpack_ILUT::Condest(const Ifpack_CondestType CT, 
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
Ifpack_ILUT::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_ILUT: " << Label() << endl << endl;
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
