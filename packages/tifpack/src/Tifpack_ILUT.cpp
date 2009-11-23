/*@HEADER
// ***********************************************************************
// 
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
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

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Preconditioner.hpp"
#include "Tifpack_ILUT.hpp"
#include "Tifpack_Condest.hpp"
#include "Tifpack_Utils.hpp"
#include "Tifpack_HashTable.hpp"
#include "Tpetra_SerialComm.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Util.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include <functional>

using namespace Teuchos;

//==============================================================================
// FIXME: allocate Comm_ and Time_ the first Initialize() call
Tifpack_ILUT::Tifpack_ILUT(const Tpetra_RowMatrix* A) :
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
Tifpack_ILUT::~Tifpack_ILUT()
{
  Destroy();
}

//==============================================================================
void Tifpack_ILUT::Destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
}

//==========================================================================
int Tifpack_ILUT::SetParameters(Teuchos::ParameterList& List)
{
  try 
  {
    LevelOfFill_ = List.get<double>("fact: ilut level-of-fill", LevelOfFill());
    if (LevelOfFill_ <= 0.0)
      TIFPACK_CHK_ERR(-2); // must be greater than 0.0

    Athresh_ = List.get<double>("fact: absolute threshold", Athresh_);
    Rthresh_ = List.get<double>("fact: relative threshold", Rthresh_);
    Relax_ = List.get<double>("fact: relax value", Relax_);
    DropTolerance_ = List.get<double>("fact: drop tolerance", DropTolerance_);

    Label_ = "TIFPACK ILUT (fill=" + Tifpack_toString(LevelOfFill())
      + ", relax=" + Tifpack_toString(RelaxValue())
      + ", athr=" + Tifpack_toString(AbsoluteThreshold())
      + ", rthr=" + Tifpack_toString(RelativeThreshold())
      + ", droptol=" + Tifpack_toString(DropTolerance())
      + ")";
    return(0);
  }
  catch (...)
  {
    cerr << "Caught an exception while parsing the parameter list" << endl;
    cerr << "This typically means that a parameter was set with the" << endl;
    cerr << "wrong type (for example, int instead of double). " << endl;
    cerr << "please check the documentation for the type required by each parameer." << endl;
    TIFPACK_CHK_ERR(-1);
  }

  return(0);
}

//==========================================================================
int Tifpack_ILUT::Initialize()
{
  // delete previously allocated factorization
  Destroy();

  Time_.ResetStartTime();

  // check only in serial
  if (Comm().NumProc() == 1 && Matrix().NumMyRows() != Matrix().NumMyCols())
    TIFPACK_CHK_ERR(-2);
    
  NumMyRows_ = Matrix().NumMyRows();

  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_.ElapsedTime();

  return(0);
}

//==========================================================================
class Tifpack_AbsComp 
{
 public:
  inline bool operator()(const double& x, const double& y) 
  {
    return(std::abs(x) > std::abs(y));
  }
};

//==========================================================================
// MS // completely rewritten the algorithm on 25-Jan-05, using more STL
// MS // and in a nicer and cleaner way. Also, it is more efficient.
// MS // WARNING: Still not fully tested!
int Tifpack_ILUT::Compute() 
{
  if (!IsInitialized()) 
    TIFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  NumMyRows_ = A_.NumMyRows();
  int Length = A_.MaxNumEntries();
  vector<int>    RowIndicesL(Length);
  vector<double> RowValuesL(Length);
  vector<int>    RowIndicesU(Length);
  vector<double> RowValuesU(Length);
  bool distributed = (Comm().NumProc() > 1)?true:false;

  if (distributed)
  {
    SerialComm_ = rcp(new Tpetra_SerialComm);
    SerialMap_ = rcp(new Tpetra_Map(NumMyRows_, 0, *SerialComm_));
    assert (SerialComm_.get() != 0);
    assert (SerialMap_.get() != 0);
  }
  else
    SerialMap_ = rcp(const_cast<Tpetra_Map*>(&A_.RowMatrixRowMap()), false);
  
  int RowNnzU;

  L_ = rcp(new Tpetra_CrsMatrix(Copy, *SerialMap_, 0));
  U_ = rcp(new Tpetra_CrsMatrix(Copy, *SerialMap_, 0));

  if ((L_.get() == 0) || (U_.get() == 0))
    TIFPACK_CHK_ERR(-5); // memory allocation error

  // insert first row in U_ and L_
  TIFPACK_CHK_ERR(A_.ExtractMyRowCopy(0,Length,RowNnzU,
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

  // modify diagonal
  for (int i = 0 ;i < RowNnzU ; ++i) {
    if (RowIndicesU[i] == 0) {
      double& v = RowValuesU[i];
      v = AbsoluteThreshold() * EPETRA_SGN(v) + RelativeThreshold() * v;
      break;
    }
  }
  
  TIFPACK_CHK_ERR(U_->InsertGlobalValues(0,RowNnzU,&(RowValuesU[0]),
                                        &(RowIndicesU[0])));
   // FIXME: DOES IT WORK IN PARALLEL ??
  RowValuesU[0] = 1.0;
  RowIndicesU[0] = 0;
  TIFPACK_CHK_ERR(L_->InsertGlobalValues(0,1,&(RowValuesU[0]),
                                        &(RowIndicesU[0])));

  int hash_size = 128;
  while (hash_size < (int) 1.5 * A_.MaxNumEntries() * LevelOfFill())
    hash_size *= 2;

  Tifpack_HashTable SingleRowU(hash_size - 1, 1);
  Tifpack_HashTable SingleRowL(hash_size - 1, 1);

  vector<int> keys;      keys.reserve(hash_size * 10);
  vector<double> values; values.reserve(hash_size * 10);
  vector<double> AbsRow; AbsRow.reserve(hash_size * 10);

  // =================== //
  // start factorization //
  // =================== //
  
  double this_proc_flops = 0.0;

  for (int row_i = 1 ; row_i < NumMyRows_ ; ++row_i) 
  {
    // get row `row_i' of the matrix, store in U pointers
    TIFPACK_CHK_ERR(A_.ExtractMyRowCopy(row_i,Length,RowNnzU,
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

    int flops = 0;
  
    for (int col_k = start_col ; col_k < row_i ; ++col_k) {

      int*    ColIndicesK;
      double* ColValuesK;
      int     ColNnzK;

      TIFPACK_CHK_ERR(U_->ExtractGlobalRowView(col_k, ColNnzK, ColValuesK, 
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
      
      double xxx = SingleRowU.get(col_k);
      if (std::abs(xxx) > DropTolerance()) {
        SingleRowL.set(col_k, xxx / DiagonalValueK);
        ++flops;

        for (int j = 0 ; j < ColNnzK ; ++j) {
          int col_j = ColIndicesK[j];

          if (col_j < col_k) continue;

          double yyy = SingleRowL.get(col_k);
          if (yyy !=  0.0)
            SingleRowU.set(col_j, -yyy * ColValuesK[j], true);
          flops += 2;
        }
      }
    }

    this_proc_flops += 1.0 * flops;

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
      if (std::abs(values[i]) > DropTolerance()) {
        AbsRow[count++] = std::abs(values[i]);
      }

    if (count > FillL) {
      nth_element(AbsRow.begin(), AbsRow.begin() + FillL, AbsRow.begin() + count, 
                  greater<double>());
      cutoff = AbsRow[FillL];
    }

    for (int i = 0; i < sizeL; ++i) {
      if (std::abs(values[i]) >= cutoff) {
        TIFPACK_CHK_ERR(L_->InsertGlobalValues(row_i,1, &values[i], (int*)&keys[i]));
      }
      else
        DiscardedElements += values[i];
    }

    // FIXME: DOES IT WORK IN PARALLEL ???
    // add 1 to the diagonal
    double dtmp = 1.0;
    TIFPACK_CHK_ERR(L_->InsertGlobalValues(row_i,1, &dtmp, &row_i));

    // same business with U_
    count = 0;
    int sizeU = SingleRowU.getNumEntries();
    AbsRow.resize(sizeU + 1);
    keys.resize(sizeU + 1);
    values.resize(sizeU + 1);

    SingleRowU.arrayify(&keys[0], &values[0]);

    for (int i = 0; i < sizeU; ++i)
      if (keys[i] >= row_i && std::abs(values[i]) > DropTolerance())
      {
        AbsRow[count++] = std::abs(values[i]);
      }

    if (count > FillU) {
      nth_element(AbsRow.begin(), AbsRow.begin() + FillU, AbsRow.begin() + count, 
                  greater<double>());
      cutoff = AbsRow[FillU];
    }

    // sets the factors in U_
    for (int i = 0; i < sizeU; ++i) 
    {
      int col = keys[i];
      double val = values[i];

      if (col >= row_i) {
        if (std::abs(val) >= cutoff || row_i == col) {
          TIFPACK_CHK_ERR(U_->InsertGlobalValues(row_i,1, &val, &col));
        }
        else
          DiscardedElements += val;
      }
    }

    // FIXME: not so sure of that!
    if (RelaxValue() != 0.0) {
      DiscardedElements *= RelaxValue();
      TIFPACK_CHK_ERR(U_->InsertGlobalValues(row_i,1, &DiscardedElements,
                                            &row_i));
    }
  }

  double tf;
  Comm().SumAll(&this_proc_flops, &tf, 1);
  ComputeFlops_ += tf;

  TIFPACK_CHK_ERR(L_->FillComplete());
  TIFPACK_CHK_ERR(U_->FillComplete());

#if 0
  // to check the complete factorization
  Tpetra_Vector LHS(A_.RowMatrixRowMap());
  Tpetra_Vector RHS1(A_.RowMatrixRowMap());
  Tpetra_Vector RHS2(A_.RowMatrixRowMap());
  Tpetra_Vector RHS3(A_.RowMatrixRowMap());
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

  int MyNonzeros = L_->NumGlobalNonzeros() + U_->NumGlobalNonzeros();
  Comm().SumAll(&MyNonzeros, &GlobalNonzeros_, 1);

  IsComputed_ = true;

  ++NumCompute_;
  ComputeTime_ += Time_.ElapsedTime();

  return(0);

}
  
//=============================================================================
int Tifpack_ILUT::ApplyInverse(const Tpetra_MultiVector& X, 
			     Tpetra_MultiVector& Y) const
{
  if (!IsComputed())
    TIFPACK_CHK_ERR(-2); // compute preconditioner first

  if (X.NumVectors() != Y.NumVectors()) 
    TIFPACK_CHK_ERR(-3); // Return error: X and Y not the same size

  Time_.ResetStartTime();

  // NOTE: L_ and U_ are based on SerialMap_, while Xcopy is based
  // on A.Map()... which are in general different. However, Solve()
  // does not seem to care... which is fine with me.
  //
  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr<const Tpetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Tpetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  if (!UseTranspose_)
  {
    // solves LU Y = X 
    TIFPACK_CHK_ERR(L_->Solve(false,false,false,*Xcopy,Y));
    TIFPACK_CHK_ERR(U_->Solve(true,false,false,Y,Y));
  }
  else
  {
    // solves U(trans) L(trans) Y = X
    TIFPACK_CHK_ERR(U_->Solve(true,true,false,*Xcopy,Y));
    TIFPACK_CHK_ERR(L_->Solve(false,true,false,Y,Y));
  }

  ++NumApplyInverse_;
  ApplyInverseFlops_ += X.NumVectors() * 2 * GlobalNonzeros_;
  ApplyInverseTime_ += Time_.ElapsedTime();

  return(0);

}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Tifpack_ILUT::Apply(const Tpetra_MultiVector& X, 
		      Tpetra_MultiVector& Y) const 
{
  return(-98);
}

//=============================================================================
double Tifpack_ILUT::Condest(const Tifpack_CondestType CT, 
                            const int MaxIters, const double Tol,
			    Tpetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // NOTE: this is computing the *local* condest
  if (Condest_ == -1.0)
    Condest_ = Tifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//=============================================================================
std::ostream&
Tifpack_ILUT::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Tifpack_ILUT: " << Label() << endl << endl;
    os << "Level-of-fill      = " << LevelOfFill() << endl;
    os << "Absolute threshold = " << AbsoluteThreshold() << endl;
    os << "Relative threshold = " << RelativeThreshold() << endl;
    os << "Relax value        = " << RelaxValue() << endl;
    os << "Condition number estimate       = " << Condest() << endl;
    os << "Global number of rows           = " << A_.NumGlobalRows() << endl;
    if (IsComputed_) {
      os << "Number of nonzeros in A         = " << A_.NumGlobalNonzeros() << endl;
      os << "Number of nonzeros in L + U     = " << NumGlobalNonzeros() 
         << " ( = " << 100.0 * NumGlobalNonzeros() / A_.NumGlobalNonzeros() 
         << " % of A)" << endl;
      os << "nonzeros / rows                 = " 
        << 1.0 * NumGlobalNonzeros() / U_->NumGlobalRows() << endl;
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
