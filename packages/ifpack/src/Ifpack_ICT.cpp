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
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Preconditioner.h"
#include "Ifpack_ICT.h"
#include "Ifpack_Condest.h"
#include "Ifpack_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Util.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
// FIXME: allocate Comm_ and Time_ the first Initialize() call
Ifpack_ICT::Ifpack_ICT(const Epetra_RowMatrix* A) :
  A_(*A),
  Comm_(Comm()),
  H_(0),
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
  Time_(Comm())
{
}

//==============================================================================
Ifpack_ICT::~Ifpack_ICT()
{
  Destroy();
}

//==============================================================================
void Ifpack_ICT::Destroy()
{
  if (H_) delete H_; H_ = 0;

  IsInitialized_ = false;
  IsComputed_ = false;
}

//==========================================================================
int Ifpack_ICT::SetParameters(Teuchos::ParameterList& List)
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
    + ")";

  return(0);
}

//==========================================================================
int Ifpack_ICT::Initialize()
{
  // clean data if present
  Destroy();

  Time_.ResetStartTime();

  // matrix must be square
  if (Matrix().NumMyRows() != Matrix().NumMyCols())
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

  int NumMyRows_ = A_.NumMyRows();
  int Length = A_.MaxNumEntries();
  vector<int>    RowIndices(Length);
  vector<double> RowValues(Length);

  int RowNnz;
  double flops = 0.0;

  H_ = new Epetra_CrsMatrix(Copy,A_.RowMatrixRowMap(),0);
  if (H_ == 0)
    IFPACK_CHK_ERR(-5);

  // get A(0,0) element and insert it (after sqrt)
  IFPACK_CHK_ERR(A_.ExtractMyRowCopy(0,Length,RowNnz,
                                     &RowValues[0],&RowIndices[0]));

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

  // start factorization for line 1
  for (int row_i = 1 ; row_i < NumMyRows_ ; ++row_i) {

    // get row `row_i' of the matrix
    IFPACK_CHK_ERR(A_.ExtractMyRowCopy(row_i,Length,RowNnz,
                                       &RowValues[0],&RowIndices[0]));

    // number of nonzeros in this row are defined as the nonzeros
    // of the matrix, plus the level of fill 
    int LOF = (int)(LevelOfFill() * RowNnz);
    if (LOF == 0) LOF = 1;

    // convert line `row_i' into STL map for fast access
    map<int,double> SingleRow;
    map<int,double>::iterator where;

    double h_ii = 0.0;
    for (int i = 0 ; i < RowNnz ; ++i) {
      if (RowIndices[i] == row_i) {
        double& v = RowValues[i];
        h_ii = AbsoluteThreshold() * EPETRA_SGN(v) + RelativeThreshold() * v;
      }
      else if (RowIndices[i] < row_i)
        SingleRow[RowIndices[i]] = RowValues[i];
    }
      
    // form element (row_i, col_j)
    // I start from the first row that has a nonzero column
    // index in row_i.
    for (int col_j = RowIndices[0] ; col_j < row_i ; ++col_j) {

      short int flops = 0;
      double h_ij = 0.0, h_jj = 0.0;
      where = SingleRow.find(col_j);
      if (where != SingleRow.end())
        h_ij = SingleRow[(*where).first];

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
          where = SingleRow.find(col_k);
          if (where != SingleRow.end()) {
            h_ij -= ColValues[k] * SingleRow[(*where).first];
            flops += 2;
          }
        }
      }

      h_ij /= h_jj;

      if (IFPACK_ABS(h_ij) > DropTolerance_)
        SingleRow[col_j] = h_ij;
    
      // only approx
      ComputeFlops_ += 2.0 * flops + 1;
    }

    int size = SingleRow.size();

    vector<double> AbsRow(size);
    int count = 0;
    for (where = SingleRow.begin() ; where != SingleRow.end() ; ++where) {
      AbsRow[count++] = IFPACK_ABS((*where).second);
    }

    double cutoff = 0.0;
    if (count > LOF) {
      nth_element(AbsRow.begin(), AbsRow.end() + LOF, AbsRow.end(), 
                  greater<double>());
      cutoff = AbsRow[LOF];
    }

    // fix the diagonal element
    for (where = SingleRow.begin() ; where != SingleRow.end() ; ++where) {
      h_ii -= ((*where).second) * ((*where).second);
    }
    if (h_ii < 0.0) h_ii = 1e-12;;

    SingleRow[row_i] = sqrt(h_ii);

    // only approx, + 1 == sqrt
    ComputeFlops_ += 2 * SingleRow.size() + 1;

    double DiscardedElements = 0.0;

    for (where = SingleRow.begin() ; where != SingleRow.end() ; ++where) {
      if (IFPACK_ABS((*where).second) > cutoff) {
        IFPACK_CHK_ERR(H_->InsertGlobalValues(row_i,1, &((*where).second),
                                              (int*)&((*where).first)));
      }
      else
        DiscardedElements += (*where).second;
    }

    // FIXME: not so sure of that!
    if (RelaxValue() != 0.0) {
      DiscardedElements *= RelaxValue();
      IFPACK_CHK_ERR(H_->InsertGlobalValues(row_i,1, &DiscardedElements,
                                            &row_i));
    }

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
#endif

  IsComputed_ = true;
  double TotalFlops; // sum across all the processors
  A_.Comm().SumAll(&flops, &TotalFlops, 1);
  ComputeFlops_ += TotalFlops;
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
  const Epetra_MultiVector* Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = new Epetra_MultiVector(X);
  else
    Xcopy = &X;

  EPETRA_CHK_ERR(H_->Solve(false,false,false,*Xcopy,Y));
  EPETRA_CHK_ERR(H_->Solve(false,true,false,Y,Y));

  if (Xcopy != &X)
    delete Xcopy;

  // these are global flop count
  ApplyInverseFlops_ += 4.0 * H_->NumGlobalNonzeros();

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
#endif
