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
#include "Ifpack_ICT.h"
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

//==============================================================================
// FIXME: allocate Comm_ and Time_ the first Initialize() call
Ifpack_ICT::Ifpack_ICT(const Epetra_RowMatrix* A) :
  A_(*A),
  Comm_(A_.Comm()),
  Condest_(-1.0),
  Athresh_(0.0),
  Rthresh_(1.0),
  LevelOfFill_(1.0),
  DropTolerance_(0.0),
  Relax_(0.0),
  IsInitialized_(false),
  IsComputed_(false),
  UseTranspose_(false),
  NumMyRows_(0),
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
  // do nothing here
}

//==============================================================================
Ifpack_ICT::~Ifpack_ICT()
{
  Destroy();
}

//==============================================================================
void Ifpack_ICT::Destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
}

//==========================================================================
int Ifpack_ICT::SetParameters(Teuchos::ParameterList& List)
{

  try
  {
    LevelOfFill_ = List.get("fact: ict level-of-fill",LevelOfFill_);
    Athresh_ = List.get("fact: absolute threshold", Athresh_);
    Rthresh_ = List.get("fact: relative threshold", Rthresh_);
    Relax_ = List.get("fact: relax value", Relax_);
    DropTolerance_ = List.get("fact: drop tolerance", DropTolerance_);

    // set label
    Label_ = "ICT (fill=" + Ifpack_toString(LevelOfFill())
      + ", athr=" + Ifpack_toString(AbsoluteThreshold()) 
      + ", rthr=" + Ifpack_toString(RelativeThreshold())
      + ", relax=" + Ifpack_toString(RelaxValue())
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
}

//==========================================================================
int Ifpack_ICT::Initialize()
{
  // clean data if present
  Destroy();

  Time_.ResetStartTime();

  // matrix must be square. Check only on one processor
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
int Ifpack_ICT::Compute() 
{
  if (!IsInitialized()) 
    IFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  NumMyRows_ = A_.NumMyRows();
  int Length = A_.MaxNumEntries();
  vector<int>    RowIndices(Length);
  vector<double> RowValues(Length);

  bool distributed = (Comm().NumProc() > 1)?true:false;

  if (distributed)
  {
    SerialComm_ = Teuchos::rcp(new Epetra_SerialComm);
    SerialMap_ = Teuchos::rcp(new Epetra_Map(NumMyRows_, 0, *SerialComm_));
    assert (SerialComm_.get() != 0);
    assert (SerialMap_.get() != 0);
  }
  else
    SerialMap_ = Teuchos::rcp(const_cast<Epetra_Map*>(&A_.RowMatrixRowMap()), false);

  int RowNnz;
#ifdef IFPACK_FLOPCOUNTERS
  double flops = 0.0;
#endif

  H_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*SerialMap_,0));
  if (H_.get() == 0)
    IFPACK_CHK_ERR(-5); // memory allocation error

  // get A(0,0) element and insert it (after sqrt)
  IFPACK_CHK_ERR(A_.ExtractMyRowCopy(0,Length,RowNnz,
                                     &RowValues[0],&RowIndices[0]));

  // skip off-processor elements
  if (distributed)
  {
    int count = 0;
    for (int i = 0 ;i < RowNnz ; ++i) 
    {
      if (RowIndices[i] < NumMyRows_){
        RowIndices[count] = RowIndices[i];
        RowValues[count] = RowValues[i];
        ++count;
      }
      else
        continue;
    }
    RowNnz = count;
  }

  // modify diagonal
  double diag_val = 0.0;
  for (int i = 0 ;i < RowNnz ; ++i) {
    if (RowIndices[i] == 0) {
      double& v = RowValues[i];
      diag_val = AbsoluteThreshold() * EPETRA_SGN(v) +
        RelativeThreshold() * v;
      break;
    }
  }

  diag_val = sqrt(diag_val);
  int diag_idx = 0;
  EPETRA_CHK_ERR(H_->InsertGlobalValues(0,1,&diag_val, &diag_idx));

  int oldSize = RowNnz;

  // The 10 is just a small constant to limit collisons as the actual keys
  // we store are the indices and not integers
  // [0..A_.MaxNumEntries()*LevelofFill()].
  Ifpack_HashTable Hash( 10 * A_.MaxNumEntries() * LevelOfFill(), 1);

  // start factorization for line 1
  for (int row_i = 1 ; row_i < NumMyRows_ ; ++row_i) {

    // get row `row_i' of the matrix
    IFPACK_CHK_ERR(A_.ExtractMyRowCopy(row_i,Length,RowNnz,
                                       &RowValues[0],&RowIndices[0]));

    // skip off-processor elements
    if (distributed)
    {
      int count = 0;
      for (int i = 0 ;i < RowNnz ; ++i) 
      {
        if (RowIndices[i] < NumMyRows_){
          RowIndices[count] = RowIndices[i];
          RowValues[count] = RowValues[i];
          ++count;
        }
        else
          continue;
      }
      RowNnz = count;
    }

    // number of nonzeros in this row are defined as the nonzeros
    // of the matrix, plus the level of fill 
    int LOF = (int)(LevelOfFill() * RowNnz);
    if (LOF == 0) LOF = 1;

    // convert line `row_i' into hash for fast access
    Hash.reset();

    double h_ii = 0.0;
    for (int i = 0 ; i < RowNnz ; ++i) {
      if (RowIndices[i] == row_i) {
        double& v = RowValues[i];
        h_ii = AbsoluteThreshold() * EPETRA_SGN(v) + RelativeThreshold() * v;
      }
      else if (RowIndices[i] < row_i)
      {
        Hash.set(RowIndices[i], RowValues[i], true);
      }
    }
      
    // form element (row_i, col_j)
    // I start from the first row that has a nonzero column
    // index in row_i.
    for (int col_j = RowIndices[0] ; col_j < row_i ; ++col_j) {

      double h_ij = 0.0, h_jj = 0.0;
      // note: get() returns 0.0 if col_j is not found
      h_ij = Hash.get(col_j);

      // get pointers to row `col_j'
      int* ColIndices;
      double* ColValues;
      int ColNnz;
      H_->ExtractGlobalRowView(col_j, ColNnz, ColValues, ColIndices);

      for (int k = 0 ; k < ColNnz ; ++k) {
        int col_k = ColIndices[k];

        if (col_k == col_j)
          h_jj = ColValues[k];
        else {
          double xxx = Hash.get(col_k);
          if (xxx != 0.0)
          {
            h_ij -= ColValues[k] * xxx;
#ifdef IFPACK_FLOPCOUNTERS
            flops += 2.0;
#endif
          }
        }
      }

      h_ij /= h_jj;

      if (IFPACK_ABS(h_ij) > DropTolerance_)
      {
        Hash.set(col_j, h_ij);
      }
    
#ifdef IFPACK_FLOPCOUNTERS
      // only approx
      ComputeFlops_ += 2.0 * flops + 1.0;
#endif
    }

    int size = Hash.getNumEntries();

    vector<double> AbsRow(size);
    int count = 0;
    
    // +1 because I use the extra position for diagonal in insert
    vector<int> keys(size + 1);
    vector<double> values(size + 1);

    Hash.arrayify(&keys[0], &values[0]);

    for (int i = 0 ; i < size ; ++i)
    {
      AbsRow[i] = IFPACK_ABS(values[i]);
    }
    count = size;

    double cutoff = 0.0;
    if (count > LOF) {
      nth_element(AbsRow.begin(), AbsRow.begin() + LOF, AbsRow.begin() + count, 
                  greater<double>());
      cutoff = AbsRow[LOF];
    }

    for (int i = 0 ; i < size ; ++i)
    {
      h_ii -= values[i] * values[i];
    }

    if (h_ii < 0.0) h_ii = 1e-12;;

    h_ii = sqrt(h_ii);

#ifdef IFPACK_FLOPCOUNTERS
    // only approx, + 1 == sqrt
    ComputeFlops_ += 2 * size + 1;
#endif

    double DiscardedElements = 0.0;

    count = 0;
    for (int i = 0 ; i < size ; ++i)    
    { 
      if (IFPACK_ABS(values[i]) > cutoff)
      {
        values[count] = values[i];
        keys[count] = keys[i];
        ++count;
      }
      else  
        DiscardedElements += values[i];
    }

    if (RelaxValue() != 0.0) {
      DiscardedElements *= RelaxValue();
      h_ii += DiscardedElements;
    }

    values[count] = h_ii;
    keys[count] = row_i;
    ++count;

    H_->InsertGlobalValues(row_i, count, &(values[0]), (int*)&(keys[0]));

    oldSize = size;
  }

  IFPACK_CHK_ERR(H_->FillComplete());

#if 0
  // to check the complete factorization
  Epetra_Vector LHS(Matrix().RowMatrixRowMap());
  Epetra_Vector RHS1(Matrix().RowMatrixRowMap());
  Epetra_Vector RHS2(Matrix().RowMatrixRowMap());
  Epetra_Vector RHS3(Matrix().RowMatrixRowMap());
  LHS.Random();

  Matrix().Multiply(false,LHS,RHS1);
  H_->Multiply(true,LHS,RHS2);
  H_->Multiply(false,RHS2,RHS3);

  RHS1.Update(-1.0, RHS3, 1.0);
  cout << endl;
  cout << RHS1;
#endif
  int MyNonzeros = H_->NumGlobalNonzeros();
  Comm().SumAll(&MyNonzeros, &GlobalNonzeros_, 1);

  IsComputed_ = true;
#ifdef IFPACK_FLOPCOUNTERS
  double TotalFlops; // sum across all the processors
  A_.Comm().SumAll(&flops, &TotalFlops, 1);
  ComputeFlops_ += TotalFlops;
#endif
  ++NumCompute_;
  ComputeTime_ += Time_.ElapsedTime();

  return(0);

}

//=============================================================================
int Ifpack_ICT::ApplyInverse(const Epetra_MultiVector& X, 
			     Epetra_MultiVector& Y) const
{

  if (!IsComputed())
    IFPACK_CHK_ERR(-3); // compute preconditioner first

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK_CHK_ERR(-2); // Return error: X and Y not the same size

  Time_.ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr<const Epetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  // NOTE: H_ is based on SerialMap_, while Xcopy is based
  // on A.Map()... which are in general different. However, Solve()
  // does not seem to care... which is fine with me.
  //
  EPETRA_CHK_ERR(H_->Solve(false,false,false,*Xcopy,Y));
  EPETRA_CHK_ERR(H_->Solve(false,true,false,Y,Y));

#ifdef IFPACK_FLOPCOUNTERS
  // these are global flop count
  ApplyInverseFlops_ += 4.0 * GlobalNonzeros_;
#endif

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_.ElapsedTime();

  return(0);
}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_ICT::Apply(const Epetra_MultiVector& X, 
		      Epetra_MultiVector& Y) const 
{

  IFPACK_CHK_ERR(-98);
}

//=============================================================================
double Ifpack_ICT::Condest(const Ifpack_CondestType CT, 
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
Ifpack_ICT::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_ICT: " << Label() << endl << endl;
    os << "Level-of-fill      = " << LevelOfFill() << endl;
    os << "Absolute threshold = " << AbsoluteThreshold() << endl;
    os << "Relative threshold = " << RelativeThreshold() << endl;
    os << "Relax value        = " << RelaxValue() << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows            = " << Matrix().NumGlobalRows() << endl;
    if (IsComputed_) {
      os << "Number of nonzeros of H         = " << H_->NumGlobalNonzeros() << endl;
      os << "nonzeros / rows                 = " 
         << 1.0 * H_->NumGlobalNonzeros() / H_->NumGlobalRows() << endl;
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
