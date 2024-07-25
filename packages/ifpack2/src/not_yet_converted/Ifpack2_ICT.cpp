// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_ICT.hpp"
#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Utils.hpp"
#include "Ifpack2_HashTable.hpp"
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

//==============================================================================
// FIXME: allocate Comm_ and Time_ the first Initialize() call
Ifpack2_ICT::Ifpack2_ICT(const Tpetra_RowMatrix* A) :
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
Ifpack2_ICT::~Ifpack2_ICT()
{
  Destroy();
}

//==============================================================================
void Ifpack2_ICT::Destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
}

//==========================================================================
int Ifpack2_ICT::SetParameters(Teuchos::ParameterList& List)
{

  try
  {
    LevelOfFill_ = List.get("fact: ict level-of-fill",LevelOfFill_);
    Athresh_ = List.get("fact: absolute threshold", Athresh_);
    Rthresh_ = List.get("fact: relative threshold", Rthresh_);
    Relax_ = List.get("fact: relax value", Relax_);
    DropTolerance_ = List.get("fact: drop tolerance", DropTolerance_);

    // set label
    Label_ = "ICT (fill=" + Ifpack2_toString(LevelOfFill())
      + ", athr=" + Ifpack2_toString(AbsoluteThreshold()) 
      + ", rthr=" + Ifpack2_toString(RelativeThreshold())
      + ", relax=" + Ifpack2_toString(RelaxValue())
      + ", droptol=" + Ifpack2_toString(DropTolerance())
      + ")";

    return(0);
  }
  catch (...)
  {
    cerr << "Caught an exception while parsing the parameter list" << endl;
    cerr << "This typically means that a parameter was set with the" << endl;
    cerr << "wrong type (for example, int instead of double). " << endl;
    cerr << "please check the documentation for the type required by each parameer." << endl;
    IFPACK2_CHK_ERR(-1);
  }
}

//==========================================================================
int Ifpack2_ICT::Initialize()
{
  // clean data if present
  Destroy();

  Time_.ResetStartTime();

  // matrix must be square. Check only on one processor
  if (Comm().NumProc() == 1 && Matrix().NumMyRows() != Matrix().NumMyCols())
    IFPACK2_CHK_ERR(-2);
    
  NumMyRows_ = Matrix().NumMyRows();

  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_.ElapsedTime();

  return(0);
}

//==========================================================================
int Ifpack2_ICT::Compute() 
{
  if (!IsInitialized()) 
    IFPACK2_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  NumMyRows_ = A_.NumMyRows();
  int Length = A_.MaxNumEntries();
  vector<int>    RowIndices(Length);
  vector<double> RowValues(Length);

  bool distributed = (Comm().NumProc() > 1)?true:false;

  if (distributed)
  {
    SerialComm_ = Teuchos::rcp(new Tpetra_SerialComm);
    SerialMap_ = Teuchos::rcp(new Tpetra_Map(NumMyRows_, 0, *SerialComm_));
    assert (SerialComm_.get() != 0);
    assert (SerialMap_.get() != 0);
  }
  else
    SerialMap_ = Teuchos::rcp(const_cast<Tpetra_Map*>(&A_.RowMatrixRowMap()), false);

  int RowNnz;
  double flops = 0.0;

  H_ = Teuchos::rcp(new Tpetra_CrsMatrix(Copy,*SerialMap_,0));
  if (H_.get() == 0)
    IFPACK2_CHK_ERR(-5); // memory allocation error

  // get A(0,0) element and insert it (after sqrt)
  IFPACK2_CHK_ERR(A_.ExtractMyRowCopy(0,Length,RowNnz,
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

  // this is a bit magic
  int hash_size = 1024;
  while (hash_size < (int)(A_.MaxNumEntries() * LevelOfFill()))
    hash_size *= 2;
  Ifpack2_HashTable Hash(hash_size - 1, 1);

  // start factorization for line 1
  for (int row_i = 1 ; row_i < NumMyRows_ ; ++row_i) {

    // get row `row_i' of the matrix
    IFPACK2_CHK_ERR(A_.ExtractMyRowCopy(row_i,Length,RowNnz,
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
            flops += 2.0;
          }
        }
      }

      h_ij /= h_jj;

      if (std::abs(h_ij) > DropTolerance_)
      {
        Hash.set(col_j, h_ij);
      }
    
      // only approx
      ComputeFlops_ += 2.0 * flops + 1.0;
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
      AbsRow[i] = std::abs(values[i]);
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

    // only approx, + 1 == sqrt
    ComputeFlops_ += 2 * size + 1;

    double DiscardedElements = 0.0;

    count = 0;
    for (int i = 0 ; i < size ; ++i)    
    { 
      if (std::abs(values[i]) > cutoff)
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

  IFPACK2_CHK_ERR(H_->FillComplete());

#if 0
  // to check the complete factorization
  Tpetra_Vector LHS(Matrix().RowMatrixRowMap());
  Tpetra_Vector RHS1(Matrix().RowMatrixRowMap());
  Tpetra_Vector RHS2(Matrix().RowMatrixRowMap());
  Tpetra_Vector RHS3(Matrix().RowMatrixRowMap());
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
  double TotalFlops; // sum across all the processors
  A_.Comm().SumAll(&flops, &TotalFlops, 1);
  ComputeFlops_ += TotalFlops;
  ++NumCompute_;
  ComputeTime_ += Time_.ElapsedTime();

  return(0);

}

//=============================================================================
int Ifpack2_ICT::ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			     Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{

  if (!IsComputed())
    IFPACK2_CHK_ERR(-3); // compute preconditioner first

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK2_CHK_ERR(-2); // Return error: X and Y not the same size

  Time_.ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RCP<const Tpetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Tpetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  // NOTE: H_ is based on SerialMap_, while Xcopy is based
  // on A.Map()... which are in general different. However, Solve()
  // does not seem to care... which is fine with me.
  //
  EPETRA_CHK_ERR(H_->Solve(false,false,false,*Xcopy,Y));
  EPETRA_CHK_ERR(H_->Solve(false,true,false,Y,Y));

  // these are global flop count
  ApplyInverseFlops_ += 4.0 * GlobalNonzeros_;

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_.ElapsedTime();

  return(0);
}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack2_ICT::Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
		      Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const 
{

  IFPACK2_CHK_ERR(-98);
}

//=============================================================================
double Ifpack2_ICT::Condest(const Ifpack2_CondestType CT, 
                            const int MaxIters, const double Tol,
			    Tpetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // NOTE: this is computing the *local* condest
  if (Condest_ == -1.0)
    Condest_ = Ifpack2_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//=============================================================================
std::ostream&
Ifpack2_ICT::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack2_ICT: " << Label() << endl << endl;
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
