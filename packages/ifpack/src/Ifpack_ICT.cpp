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
  Relax_(1.0),
  Athresh_(0.0),
  Rthresh_(0.0),
  LevelOfFill_(1.0),
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
  Time_(Comm())
{
}

//==============================================================================
Ifpack_ICT::~Ifpack_ICT()
{

  if (H_)
    delete H_;

  IsInitialized_ = false;
  IsComputed_ = false;
}

//==========================================================================
int Ifpack_ICT::SetParameters(Teuchos::ParameterList& List)
{

  LevelOfFill_ = List.get("fact: ict level-of-fill",LevelOfFill_);
  Athresh_ = List.get("fact: absolute threshold", Athresh_);
  Rthresh_ = List.get("fact: relative threshold", Rthresh_);

  // set label
  sprintf(Label_, "ICT (fill=%f, athr=%f, rthr=%f)",
	  LevelOfFill_, Athresh_, 
	  Rthresh_);
  return(0);
}

//==========================================================================
int Ifpack_ICT::Initialize()
{
  IsInitialized_ = false;
  Time_.ResetStartTime();

  if (Matrix().NumMyRows() != Matrix().NumMyCols())
    IFPACK_CHK_ERR(-2);
    
  NumMyRows_ = Matrix().NumMyRows();

  // delete previously allocated factorization
  if (H_)
    delete H_;

  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_.ElapsedTime();

  return(0);
}

//==========================================================================
int Ifpack_ICT::Compute() {

  if (!IsInitialized()) 
    IFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  // extract information from matrix A_, put them
  // in STL vectors that will be used in ComputeZeroFill()
  // and ComputeNonZeroFill()

  int NumMyRows_ = Matrix().NumMyRows();

  vector<int> Nnz(NumMyRows_);
  vector<vector<int> > Indices(NumMyRows_);
  vector<vector<double> > Values(NumMyRows_);

  int Length = Matrix().MaxNumEntries();
  vector<int> RowIndices(Length);
  vector<double> RowValues(Length);

  vector<int> Diag(NumMyRows_);

  long int flops = 0;

  // first cycle over all rows and extract Nnz, indices and values
  // for each row
  //
  for (int i = 0 ; i < NumMyRows_ ; ++i) {

    int RowNnz;
    Matrix().ExtractMyRowCopy(i,Length,RowNnz,&RowValues[0],&RowIndices[0]);
    
#if 0
    for (int j = 0 ; j < RowNnz ; ++j) {
      cout << "A(" <<  1 +i << ", " << 1 + RowIndices[j] << ") = " <<
	RowValues[j] << ";" << endl;
    }
    cout << endl;
#endif

    // count how many elements are below diagonal, and store
    // the position of diagonal element
    int count = 0;
    double DiagVal = 0.0;
    for (int j = 0 ; j < RowNnz ; ++j) {
      if (RowIndices[j] == i)
	DiagVal = RowValues[j];
      else if (RowIndices[j] < i)
	++count;
    }

    ++count; // diagonal value
    Nnz[i] = count;
    Indices[i].resize(count);
    Values[i].resize(count);

    // I suppose that the matrix has diagonal value, and I store
    // it last. Matrix always has the diagonal entry (albeit possibly
    // zero) if this matrix comes from Ifpack_AdditiveSchwarz.
    count = 0;
    for (int j = 0 ; j < RowNnz ; ++j) {
      if (RowIndices[j] < i) {
	Indices[i][count] = RowIndices[j]; 
	Values[i][count] = RowValues[j]; 
	++count;
      }
    }
    Indices[i][count] = i;
    Values[i][count] = DiagVal;
    Diag[i] = count;
    ++count;

    assert (count == Nnz[i]);
  }

  H_ = new Epetra_CrsMatrix(Copy,Matrix().RowMatrixRowMap(),0);
  if (H_ == 0)
    IFPACK_CHK_ERR(-5);

  // insert element Matrix()(0,0)
  Values[0][0] = sqrt(Values[0][0]);
  EPETRA_CHK_ERR(H_->InsertGlobalValues(0,Nnz[0],&(Values[0][0]),
                                        &(Indices[0][0])));

  vector<double> tmp(NumMyRows_);

  // now perform the factorization for row `row_i'
  for (int row_i = 1 ; row_i < NumMyRows_ ; ++row_i) {

    // number of nonzeros in this row are defined as the nonzeros
    // of the matrix, plus the level of fill (-1, because Nnz
    // contains the diagonal as well, which is always stored).
    int LOF = (int)(LevelOfFill_ * (Nnz[row_i] - 1));

    // first non-zero index in row_i
    int MinIndex = Indices[row_i][0];
    // get diagonal value for atresh_
    double DiagonalValue = Values[row_i][Diag[row_i]];

    // zero-out tmp from MinIndex to subdiagonal entry
    for (int j = MinIndex ; j < row_i - 1 ; ++j)
      tmp[j] = 0.0;

    // put nonzero values of row_i
    for (int j = 0 ; j < Nnz[row_i] - 1 ; ++j) {
      int col_j = Indices[row_i][j];
      tmp[col_j] = Values[row_i][j];
    }

    // put all elements in this vector, because it
    // will be sorted using STL sort() algorithm later
    vector<Ifpack_Element> RowEntries;
    Ifpack_Element Element;

    // form element (row_i, col_j) -- all of them.
    // I start from the first row that has a nonzero column
    // index in row_i.
    for (int col_j = MinIndex ; col_j < row_i ; ++col_j) {
      double& h_ij = tmp[col_j];

      // get pointers to row `col_j'
      int* ColIndices;
      double* ColValues;
      int ColNnz;
      H_->ExtractGlobalRowView(col_j, ColNnz, ColValues, ColIndices);
      // all cross-products will be zero, skip all
      if (ColIndices[ColNnz - 1] < MinIndex)
        continue;

      // look for cross product between `row_i' and `col_j'
      for (int k = 0 ; k < ColNnz - 1 ; ++k) {
        int col_k = ColIndices[k];
        if (col_k < MinIndex)
          continue;
        h_ij -= tmp[col_k] * ColValues[k];
        flops += 2;
      }

      if (h_ij == 0.0 || IFPACK_ABS(h_ij) < Athresh_
          || IFPACK_ABS(h_ij < DiagonalValue) < Rthresh_)
        continue;

      h_ij /= ColValues[ColNnz - 1];
      ++flops;
      Element.SetIndex(col_j);
      Element.SetValue(h_ij);
      RowEntries.push_back(Element);
    }

    if ((int)RowEntries.size() > LOF) {
      // sort in ascending order the entries for this row
      sort(RowEntries.begin(),RowEntries.end());
    }

    // look for the components that are in the level-of-fill range
    // NOTE: here level-of-fill refers to the number of
    // *off-diagonal* elements. If zero, no off-diagonal elements
    // will be mantained.
    for (int k = 0 ; k < EPETRA_MIN(LOF,(int)RowEntries.size()) ; ++k) {
      int col_k = RowEntries[k].Index();
      if (col_k >= row_i)
        continue;
      double h_ij = RowEntries[k].Value();
      // skip zero elements
      if (h_ij == 0.0)
        continue;
      DiagonalValue -= h_ij * h_ij;
      flops += 2;
      H_->InsertGlobalValues(row_i,1, &h_ij, &col_k);
    }
    // FIXME: use relax here??

    // diagonal element is always inserted
    DiagonalValue = sqrt(DiagonalValue);
    H_->InsertGlobalValues(row_i,1, &DiagonalValue, &row_i);

    Values[row_i].resize(0);
    Indices[row_i].resize(0);
  }

  IFPACK_CHK_ERR(H_->FillComplete());

  IsComputed_ = true;
  ComputeFlops_ += flops;
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

  // NOTE: I do suppose that X and Y are two different vectors
  EPETRA_CHK_ERR(H_->Solve(false,false,false,X,Y));
  EPETRA_CHK_ERR(H_->Solve(false,true,false,Y,Y));

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
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops() 
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops() / ComputeTime() << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse() 
       << "  " << std::setw(15) << ApplyInverseTime()
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops() 
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops() / ApplyInverseTime() << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  
  return(os);
}
#endif
