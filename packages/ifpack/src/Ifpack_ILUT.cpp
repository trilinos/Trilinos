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
#include "Ifpack_ILUT.h"
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
#include <functional>

//==============================================================================
// FIXME: allocate Comm_ and Time_ the first Initialize() call
Ifpack_ILUT::Ifpack_ILUT(const Epetra_RowMatrix* A) :
  A_(*A),
  Comm_(Comm()),
  L_(0),
  U_(0),
  Condest_(-1.0),
  Athresh_(0.0),
  Rthresh_(0.0),
  LevelOfFill_(1.0),
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
  Time_(Comm())
{
}

//==============================================================================
Ifpack_ILUT::Ifpack_ILUT(const Ifpack_ILUT& rhs) :
  A_(rhs.Matrix()),
  Comm_(Comm()),
  L_(0),
  U_(0),
  Condest_(rhs.Condest()),
  Athresh_(rhs.AbsoluteThreshold()),
  Rthresh_(rhs.RelativeThreshold()),
  LevelOfFill_(rhs.LevelOfFill()),
  IsInitialized_(rhs.IsInitialized()),
  IsComputed_(rhs.IsComputed()),
  NumInitialize_(rhs.NumInitialize()),
  NumCompute_(rhs.NumCompute()),
  NumApplyInverse_(rhs.NumApplyInverse()),
  InitializeTime_(rhs.InitializeTime()),
  ComputeTime_(rhs.ComputeTime()),
  ApplyInverseTime_(rhs.ApplyInverseTime()),
  Time_(Comm())
{
  L_ = new Epetra_CrsMatrix(rhs.L());
  U_ = new Epetra_CrsMatrix(rhs.U());
}

//==============================================================================
Ifpack_ILUT::~Ifpack_ILUT()
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
int Ifpack_ILUT::SetParameters(Teuchos::ParameterList& List)
{

  LevelOfFill_ = List.get("fact: ilut level-of-fill",LevelOfFill());
  if (LevelOfFill_ <= 0.0)
    IFPACK_CHK_ERR(-1); // must be greater than 0.0

  Athresh_ = List.get("fact: absolute threshold", Athresh_);
  Rthresh_ = List.get("fact: relative threshold", Rthresh_);

  // set label
  sprintf(Label_, "IFPACK ILUT (fill=%f, athr=%f, rthr=%f)",
	  LevelOfFill(), AbsoluteThreshold(), 
	  RelativeThreshold());
  return(0);
}
#endif

//==========================================================================
int Ifpack_ILUT::Initialize()
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
#include "float.h"
int Ifpack_ILUT::Compute() {

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

  // sizes of L and U pointers (a couple for each row)
  vector<int> L_size(NumMyRows_);
  vector<int> U_size(NumMyRows_);

  // cycle over all rows
  for (int i = 0 ; i < NumMyRows_ ; ++i) {

    int RowNnz;
    Matrix().ExtractMyRowCopy(i,Length,RowNnz,&RowValues[0],&RowIndices[0]);
    
    // count the nonzeros in L (without diagonal, stored in U)
    L_size[i] = 0;
    U_size[i] = 0;
    for (int j = 0 ; j < RowNnz ; ++j) {
      if (RowIndices[j] < i)
        L_size[i]++;
      else
        U_size[i]++;
    }

    LI[i].resize((int)(LevelOfFill() * L_size[i]));
    LV[i].resize((int)(LevelOfFill() * L_size[i]));
    UI[i].resize((int)(LevelOfFill() * (U_size[i] + 1)));
    UV[i].resize((int)(LevelOfFill() * (U_size[i] + 1)));

    int L_count = 0;
    int U_count = 0;

    // zero diagonal is always the first of U
    UI[i][U_count] = i;
    UV[i][U_count] = 0.0;
    ++U_count;

    for (int j = 0 ; j < RowNnz ; ++j) {
      int col = RowIndices[j];

      if (col < i) {
        LI[i][L_count] = col;
        LV[i][L_count] = RowValues[j];
        ++L_count;
      }
      else {
        if (col == i) {
          UV[i][0] = RowValues[j];
        }
        else {
          UI[i][U_count] = col;
          UV[i][U_count] = RowValues[j];
          ++U_count;
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
  IFPACK_CHK_ERR(U_->InsertGlobalValues(0,(int)U_size[0],
                                        &(UV[0][0]), &(UI[0][0])));

  // FIXME
  //double threshold_ = 0e-6;
  double rel_threshold_ = 1e-1;

  // ===================== //
  // Perform factorization //
  // ===================== //

  vector<double> tmp(NumMyRows_);
  for (int j = 0 ; j < NumMyRows_ ; ++j)
    tmp[j] = 0.0;
  vector<double> l_atmp(NumMyRows_);
  vector<double> u_atmp(NumMyRows_);
  vector<int> l_index(NumMyRows_);
  vector<int> u_index(NumMyRows_);

  double old_l_cutoff = 1.0;
  double old_u_cutoff = 1.0;

  // cycle over all rows 
  for (int i = 1 ; i < NumMyRows_ ; ++i) {

    // populate tmp with nonzeros of this row
    for (int j = 0 ; j < (int)L_size[i] ; ++j) {
      tmp[LI[i][j]] = LV[i][j];
    }
    for (int j = 0 ; j < (int)U_size[i] ; ++j) {
      tmp[UI[i][j]] = UV[i][j];
    }

    int first;
    if (L_size[i])
      first = LI[i][0];
    else
      first = i;

    double diag = UV[i][0];
    double drop = 1e-4 * diag;

    for (int k = first ; k < i ; ++k) {

      if (tmp[k] == 0.0)
        continue;

      if (IFPACK_ABS(tmp[k]) < drop) {
        tmp[k] = 0.0;
        continue;
      }

      tmp[k] /= UV[k][0];

      for (int j = 0 ; j < (int)U_size[k] ; ++j) {
        int col = UI[k][j];
        if (col <= k)
          continue;
        double add = tmp[k] * UV[k][j];
        if (IFPACK_ABS(add) > drop)
          tmp[col] -= add;
      }
    }

    // track diagonal element and insert it
    if (tmp[i] < 1-9)
      tmp[i] = 1.0; // FIXME
    UV[i][0] = tmp[i];
    double abs_diag = IFPACK_ABS(UV[i][0]);
    tmp[i] = 0.0;

    // estimate a good cut-off from previous line. This will
    // limitate the number of elements to order with sort().
    double this_l_cutoff = rel_threshold_ * old_l_cutoff * abs_diag;
    double this_u_cutoff = rel_threshold_ * old_u_cutoff * abs_diag;
    // get nonzeros in tmp, and absolute values in l_atmp and u_atmp
    int l_count = 0;
    int u_count = 0;
    for (int j = 0 ; j < i ; ++j) {
      double val = IFPACK_ABS(tmp[j]);
      if (val <= this_l_cutoff || val <= drop) {
        tmp[j] = 0.0;
        continue;
      }
      // store in l pointer
      l_atmp[l_count] = val;
      l_index[l_count] = j;
      ++l_count;
    }
    for (int j = i + 1 ; j < NumMyRows_ ; ++j) {
      double val = IFPACK_ABS(tmp[j]);
      if (val <= this_u_cutoff || val <= drop) {
        tmp[j] = 0.0;
        continue;
      }
      u_atmp[u_count] = val;
      u_index[u_count] = j;
      ++u_count;
    }

    int l_LOF = (int)(LevelOfFill() * L_size[i]);
    int u_LOF = (int)(LevelOfFill() * U_size[i]);
    double l_cutoff = 0.0;
    double u_cutoff = 0.0;
    // sort in ascending order the entries for this row
    if (l_count > l_LOF) {
      sort(l_atmp.begin(),l_atmp.begin() + l_count,greater<double>());
      l_cutoff = l_atmp[l_LOF];
    }
    if (u_count > u_LOF) {
      sort(u_atmp.begin(),u_atmp.begin() + u_count,greater<double>());
      u_cutoff = u_atmp[u_LOF];
    }

    int L_count = 0;
    int U_count = 1; // diagonal already inserted

    for (int kk = 0 ; kk < l_count ; ++kk) {
      int col = l_index[kk];
      double aval = IFPACK_ABS(tmp[col]);
      if (aval > l_cutoff) {
        LI[i][L_count] = col;
        LV[i][L_count] = tmp[col];
        ++L_count;
      }
      tmp[col] = 0.0;
    }
    for (int kk = 0 ; kk < u_count ; ++kk) {
      int col = u_index[kk];
      double aval = IFPACK_ABS(tmp[col]);
      if (aval > u_cutoff) {
        UI[i][U_count] = col;
        UV[i][U_count] = tmp[col];
        ++U_count;
      }
      tmp[col] = 0.0;
    }

    // reset the number in processed row
    L_size[i] = L_count;
    U_size[i] = U_count;
    if (L_size[i] > (int)LI[i].size())
      IFPACK_CHK_ERR(-1);
    if (U_size[i] > (int)UI[i].size())
      IFPACK_CHK_ERR(-1);

    old_l_cutoff = l_cutoff / abs_diag;
    old_u_cutoff = u_cutoff / abs_diag;

  }

  // insert unit diagonal in L_
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    double val = 1.0;
    IFPACK_CHK_ERR(L_->InsertGlobalValues(i,1,&val, &i));
  }

  // insert computed elements in L_
  for (int i = 1 ; i < NumMyRows_ ; ++i) {
#ifdef IFPACK_DEBUG
    for (int j = 0 ; j < (int)L_size[i] ; ++j) {
      if (LI[i][j] >= NumMyRows_) {
        cerr << "ERROR: LI[" << i << "][" << j << "] = "
          << LI[i][j] << " and NumMyRows = " << NumMyRows_ << endl;
        cerr << "(file " << __FILE__ << ", line " << __LINE__ << endl;
        exit(EXIT_FAILURE);
      }
    }
#endif
    IFPACK_CHK_ERR(L_->InsertGlobalValues(i,(int)L_size[i],
                                          &(LV[i][0]), &(LI[i][0])));
    LI[i].resize(0);
    LV[i].resize(0);
  }

  // insert computed elements in U_
  for (int i = 1 ; i < NumMyRows_ ; ++i) {
#ifdef IFPACK_DEBUG
    for (int j = 0 ; j < (int)U_size[i] ; ++j) {
      if (UI[i][j] >= NumMyRows_) {
        cerr << "ERROR: UI[" << i << "][" << j << "] = "
          << UI[i][j] << " and NumMyRows = " << NumMyRows_ << endl;
        cerr << "(file " << __FILE__ << ", line " << __LINE__ << endl;
        exit(EXIT_FAILURE);
      }
    }
#endif
    IFPACK_CHK_ERR(U_->InsertGlobalValues(i,(int)U_size[i],
                                          &(UV[i][0]), &(UI[i][0])));
    UI[i].resize(0);
    UV[i].resize(0);
  }

  IFPACK_CHK_ERR(L_->FillComplete());
  IFPACK_CHK_ERR(U_->FillComplete());

  IsComputed_ = true;
  ++NumCompute_;
  
  return(0);

}

//=============================================================================
int Ifpack_ILUT::ApplyInverse(const Epetra_MultiVector& X, 
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

  return(0);

}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_ILUT::Apply(const Epetra_MultiVector& X, 
		      Epetra_MultiVector& Y) const 
{

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK_CHK_ERR(-1); // Return error: X and Y not the same size

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
double Ifpack_ILUT::Condest(const Ifpack_CondestType CT, 
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
Ifpack_ILUT::Print(std::ostream& os) const
{
  os << endl << "*** Ifpack_ILUT : " << Label() << endl << endl;
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
