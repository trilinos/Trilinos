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

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_vRiluk.h"
#include "Ifpack_Condest.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Util.h"

#ifdef HAVE_IFPACK_TEUCHOS
#include "Teuchos_ParameterList.hpp"
#endif

//==============================================================================
class Ifpack_Element {

public:
  Ifpack_Element() {};

  Ifpack_Element(const Ifpack_Element& rhs) {
    i_ = rhs.Index();
    val_ = rhs.Value();
    aval_ = rhs.AbsValue();
  }

  inline int Index() const {
    return(i_);
  }

  inline double Value() const {
    return(val_);
  }

  inline double AbsValue() const {
    return(aval_);
  }

  inline void SetIndex(const int i)
  {
    i_ = i;
  }

  inline void SetValue(const double val)
  {
    val_ = val;
    aval_ = IFPACK_ABS(val_);
  }

  inline bool operator <(const Ifpack_Element& rhs) const 
  {
    if (rhs.AbsValue() > AbsValue())
      return(false);
    else if (rhs.AbsValue() < AbsValue())
      return(true);
    else if (rhs.Index() < Index())
        return(true);
  }

private:
  int i_;
  double val_;
  double aval_;

};

//==============================================================================
// FIXME: allocate Comm_ and Time_ the first Initialize() call
Ifpack_vRiluk::Ifpack_vRiluk(const Epetra_RowMatrix* A) :
  A_(*A),
  Comm_(Comm()),
  Time_(Comm()),
  IsInitialized_(false),
  IsComputed_(false),
  Condest_(-1.0),
  LevelOfFill_(0),
  L_(0),
  U_(0),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  Athresh_(0.0),
  Rthresh_(0.0)
{
}

//==============================================================================
Ifpack_vRiluk::Ifpack_vRiluk(const Ifpack_vRiluk& rhs) :
  A_(rhs.Matrix()),
  Comm_(Comm()),
  Time_(Comm()),
  IsInitialized_(rhs.IsInitialized()),
  IsComputed_(rhs.IsComputed()),
  Condest_(rhs.Condest()),
  LevelOfFill_(rhs.LevelOfFill()),
  L_(0),
  U_(0),
  NumInitialize_(rhs.NumInitialize()),
  NumCompute_(rhs.NumCompute()),
  NumApplyInverse_(rhs.NumApplyInverse()),
  InitializeTime_(rhs.InitializeTime()),
  ComputeTime_(rhs.ComputeTime()),
  ApplyInverseTime_(rhs.ApplyInverseTime()),
  Athresh_(rhs.AbsoluteThreshold()),
  Rthresh_(rhs.RelativeThreshold())
{
  L_ = new Epetra_CrsMatrix(rhs.L());
  U_ = new Epetra_CrsMatrix(rhs.U());
}

//==============================================================================
Ifpack_vRiluk::~Ifpack_vRiluk()
{

  if (L_)
    delete L_;

  if (U_)
    delete U_;

  IsInitialized_ = false;
  IsComputed_ = false;
}

#ifdef HAVE_IFPACK_TEUCHOS
//==========================================================================
int Ifpack_vRiluk::SetParameters(Teuchos::ParameterList& List)
{

  LevelOfFill_ = List.get("fact: level-of-fill",LevelOfFill());
  Athresh_ = List.get("fact: absolute threshold", Athresh_);
  Rthresh_ = List.get("fact: relative threshold", Rthresh_);

  // set label
  sprintf(Label_, "vRILUK (fill=%d, athr=%f, rthr=%f)",
	  LevelOfFill(), AbsoluteThreshold(), 
	  RelativeThreshold());
  return(0);
}
#endif

//==========================================================================
int Ifpack_vRiluk::Initialize()
{
  IsInitialized_ = false;
  Time_.ResetStartTime();

  // I work on one process only (use Ifpack_AdditiveSchwarz
  // if using more than one process)
  if (Matrix().Comm().NumProc() != 1)
    IFPACK_CHK_ERR(-2);

  if (Matrix().NumMyRows() != Matrix().NumMyCols())
    IFPACK_CHK_ERR(-1);
    
  NumMyRows_ = Matrix().NumMyRows();

  // delete previously allocated factorization
  if (L_)
    delete L_;
  if (U_)
    delete U_;

  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_.ElapsedTime();

  return(0);
}

//==========================================================================
int Ifpack_vRiluk::Compute() {

  if (!IsInitialized()) 
    IFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  int NumMyRows_ = Matrix().NumMyRows();

  vector<int> Nnz(NumMyRows_);
  // will contain the nonzero values of U (including diagonal)
  vector<vector<int> > UI(NumMyRows_);
  vector<vector<double> > UV(NumMyRows_);

  // will contain the nonzero values of L
  vector<vector<int> > LI(NumMyRows_);
  vector<vector<double> > LV(NumMyRows_);

  // temp vectors for ExtractMyRowCopy()
  int Length = Matrix().MaxNumEntries();
  vector<int> RowIndices(Length);
  vector<double> RowValues(Length);

  // cycle over all rows
  for (int i = 0 ; i < NumMyRows_ ; ++i) {

    int RowNnz;
    Matrix().ExtractMyRowCopy(i,Length,RowNnz,&RowValues[0],&RowIndices[0]);
    
    // zero diagonal is always the first of U
    UI[i].push_back(i);
    UV[i].push_back(0.0);

    for (int j = 0 ; j < RowNnz ; ++j) {
      int col = RowIndices[j];

      if (col < i) {
        LI[i].push_back(col);
        LV[i].push_back(RowValues[j]);
      }
      else {
        if (col == i) {
          UV[i][0] = RowValues[j];
        }
        else {
          UI[i].push_back(col);
          UV[i].push_back(RowValues[j]);
        }
      }
    }
  }

  // allocate the Crs matrices that will contain the factors
  L_ = new Epetra_CrsMatrix(Copy,Matrix().RowMatrixRowMap(),0);
  U_ = new Epetra_CrsMatrix(Copy,Matrix().RowMatrixRowMap(),0);
  if ((L_ == 0) || (U_ == 0))
    IFPACK_CHK_ERR(-1);

  // insert first row in U_
  IFPACK_CHK_ERR(U_->InsertGlobalValues(0,UI[0].size(),
                                        &(UV[0][0]), &(UI[0][0])));

  vector<int> flag(NumMyRows_);
  for (int j = 0 ; j < NumMyRows_ ; ++j)
    flag[j] = -1;

  if (!LevelOfFill_) {

    // ================= //
    // zero fill-in case //
    // ================= //

    // cycle over all rows 
    for (int i = 1 ; i < NumMyRows_ ; ++i) {

      for (int k = 0 ; k < LI[i].size() ; ++k) {

        // index of entry `k'
        int kk = LI[i][k];

        LV[i][k] /= UV[kk][0];
        double Lik = LV[i][k];

        // position of nonzeros in row `kk'
        for (int j = 0 ; j < UI[kk].size() ; ++j)
          flag[UI[kk][j]] = j;

        // fix lower triangular
        for (int j = 0 ; j < LI[i].size() ; ++j) {
          int jj = LI[i][j];
          if (jj <= kk)
            continue;
          // look if position `jj' is in row `kk'
          if (flag[jj] != -1) {
            LV[i][j] -= Lik * UV[kk][flag[jj]];
          }
        }
        // fix upper triangular
        for (int j = 0 ; j < UI[i].size() ; ++j) {
          int jj = UI[i][j];
          if (flag[jj] != -1) {
            UV[i][j] -= Lik * UV[kk][flag[jj]];
          }
        }

        for (int j = 0 ; j < UI[kk].size() ; ++j)
          flag[UI[kk][j]] = -1;
      }
    }
  }
  else {

    // ===================== //
    // non-zero fill-in case //
    // ===================== //

    vector<double> tmp(NumMyRows_);

    // cycle over all rows 
    for (int i = 1 ; i < NumMyRows_ ; ++i) {

      // FIXME: I am not technically correct, this may keep
      // LevelOfFill() elements in L and LevelOfFill() in U...
      int LOF = 2* LevelOfFill();

      // put all elements in this vector, because it
      // will be sorted using STL sort() algorithm later
      vector<Ifpack_Element> RowEntries;
      Ifpack_Element Element;

      // populate tmp with nonzeros of this row
      for (int j = 0 ; j < NumMyRows_ ; ++j)
        tmp[j] = 0.0;
      for (int j = 0 ; j < LI[i].size() ; ++j)
        tmp[LI[i][j]] = LV[i][j];
      for (int j = 0 ; j < UI[i].size() ; ++j)
        tmp[UI[i][j]] = UV[i][j];

      // now perform factorization as if it were a dense matrix
      for (int k = 0 ; k < i ; ++k) {

        if (tmp[k] == 0.0)
          continue;

        tmp[k] /= UV[k][0];

        // position of nonzeros in row `k'
        for (int j = 0 ; j < UI[k].size() ; ++j) {
          if (UI[k][j] <= k)
            continue;
          tmp[UI[k][j]] -= tmp[k] * UV[k][j];
        }
      }

      // reset the entries in LV
      for (int j = 0 ; j < LV[i].size() ; ++j) {
        int col = LI[i][j];
        LV[i][j] = tmp[col];
        tmp[col] = 0.0;
      }
      
      // reset the entries in UV
      for (int j = 0 ; j < UV[i].size() ; ++j) {
        int col = UI[i][j];
        UV[i][j] = tmp[col];
        tmp[col] = 0.0;
      }
      
      // get nonzeros in tmp
      for (int j = 0 ; j < NumMyRows_ ; ++j) {
        if ((tmp[j] != 0.0) && (j != i)) {
          Element.SetIndex(j);
          Element.SetValue(tmp[j]);
          RowEntries.push_back(Element);
        }
      }

      if (RowEntries.size() > LOF) {
        // sort in ascending order the entries for this row
        sort(RowEntries.begin(),RowEntries.end());
      }

      for (int kk = 0 ; kk < EPETRA_MIN(LOF,RowEntries.size()) ; ++kk) {
        int col = RowEntries[kk].Index();
        double val = RowEntries[kk].Value();
        assert (val != 0.0);
        if (col < i) {
          LI[i].push_back(col);
          LV[i].push_back(val);
        }
        else {
          UI[i].push_back(col);
          UV[i].push_back(val);
        }
      }
    }
  }

  // insert unit diagonal in L_
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    double val = 1.0;
    IFPACK_CHK_ERR(L_->InsertGlobalValues(i,1,&val, &i));
  }

  // insert computed elements in L_
  for (int i = 1 ; i < NumMyRows_ ; ++i) {
    IFPACK_CHK_ERR(L_->InsertGlobalValues(i,LI[i].size(),
                                          &(LV[i][0]), &(LI[i][0])));
    LI[i].resize(0);
    LV[i].resize(0);
  }

  // insert computed elements in U_
  for (int i = 1 ; i < NumMyRows_ ; ++i) {
    IFPACK_CHK_ERR(U_->InsertGlobalValues(i,UI[i].size(),
                                          &(UV[i][0]), &(UI[i][0])));
    UI[i].resize(0);
    UV[i].resize(0);
  }

  IFPACK_CHK_ERR(L_->FillComplete());
  IFPACK_CHK_ERR(U_->FillComplete());

  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += Time_.ElapsedTime();

  return(0);

}

//=============================================================================
int Ifpack_vRiluk::ApplyInverse(const Epetra_MultiVector& X, 
			     Epetra_MultiVector& Y) const
{

  Time_.ResetStartTime();

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK_CHK_ERR(-1); // Return error: X and Y not the same size

  if (!IsComputed())
    IFPACK_CHK_ERR(-1); // compute preconditioner first

  // NOTE: I do suppose that X and Y are two different vectors
  EPETRA_CHK_ERR(L_->Solve(false,false,false,X,Y));
  EPETRA_CHK_ERR(U_->Solve(true,false,false,Y,Y));

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_.ElapsedTime();

}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_vRiluk::Apply(const Epetra_MultiVector& X, 
		      Epetra_MultiVector& Y) const 
{

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK_CHK_ERR(-1); // Return error: X and Y not the same size

  bool Upper = true;
  bool Lower = false;
  bool UnitDiagonal = true;

#ifdef FIXME
  Epetra_MultiVector * X1 = (Epetra_MultiVector *) &X;
  Epetra_MultiVector * Y1 = (Epetra_MultiVector *) &Y;

  U_->Multiply(false, *X1, *Y1);
  Y1->Update(1.0, *X1, 1.0); // Y1 = Y1 + X1 (account for implicit unit diagonal)
  Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
  Epetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
  U_->Multiply(true, Y1temp, *Y1);
  Y1->Update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
  return(0);
#endif
  return(-1);
}

//=============================================================================
double Ifpack_vRiluk::Condest(const Ifpack_CondestType CT, 
                            const int MaxIters, const double Tol,
			    Epetra_RowMatrix* Matrix)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // NOTE: this is computing the *local* condest
  if (Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix);

  return(Condest_);
}

//=============================================================================
std::ostream&
Ifpack_vRiluk::Print(std::ostream& os) const
{
  os << endl << "*** Ifpack_vRiluk : " << Label() << endl << endl;
  os << "Level-of-fill      = " << LevelOfFill() << endl;
  os << "Absolute threshold = " << AbsoluteThreshold() << endl;
  os << "Relative threshold = " << RelativeThreshold() << endl;
//  os << "Relaxation value   = " << RelaxValue() << endl;
  if (IsInitialized()) {
    os << "Preconditioner has been initialized" << endl;
  }
  if (IsComputed()) {
    os << "Preconditioner has been computed" << endl;
    os << endl;
    os << "Number of rows of L + U         = " << L_->NumMyRows() << endl;
    os << "Number of nonzeros of L + U     = " << NumGlobalNonzeros() << endl;
    os << "nonzeros / rows                 = " 
       << 1.0 * NumGlobalNonzeros() / U_->NumMyRows() << endl;
    Epetra_Vector Diagonal(U_->RowMatrixRowMap());
    U_->ExtractDiagonalCopy(Diagonal);
    double MinValue;
    double MaxValue;
    Diagonal.MinValue(&MinValue);
    Diagonal.MaxValue(&MaxValue);
    os << "Minimum diagonal value          = " << MinValue << endl;
    os << "Maximum diagonal value          = " << MaxValue << endl;
  }
  os << endl;
  os << "Number of initialization phases = " << NumInitialize_ << endl;
  os << "Number of computation phases    = " << NumCompute_ << endl;
  os << "Number of applications          = " << NumApplyInverse_ << endl;
  os << endl;
  os << "Total time for Initialize()     = " << InitializeTime_ << " (s)\n";
  os << "Total time for Compute()        = " << ComputeTime_ << " (s)\n";
  os << "Total time for ApplyInverse()   = " << ApplyInverseTime_ << " (s)\n";
  os << endl;
  
  return(os);
}



/*
 *
        // position of nonzeros in row `k'
        for (int j = 0 ; j < UI[k].size() ; ++j)
          flag[UI[k][j]] = j;

        // perform factorizations (both L and U stored in tmp)
        for (int j = k + 1 ; j < NumMyRows_ ; ++j) {
          // if I have a nonzero in this position on row `k', multiply
          if (tmp[k] != 0.0 && flag[j] != -1) {
            tmp[j] -= tmp[k] * UV[k][flag[j]];
          }
        }
*/
