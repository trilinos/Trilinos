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
#include "Ifpack_vIct.h"
#include "Ifpack_Condest.h"
#include "Ifpack_Utils.h"
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
// FIXME: allocate Comm_ and Time_ the first Initialize() call
Ifpack_vIct::Ifpack_vIct(const Epetra_RowMatrix* A) :
  A_(*A),
  Comm_(Comm()),
  Time_(Comm()),
  IsInitialized_(false),
  IsComputed_(false),
  Condest_(-1.0),
  LevelOfFill_(0),
  H_(0),
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
Ifpack_vIct::Ifpack_vIct(const Ifpack_vIct& rhs) :
  A_(rhs.Matrix()),
  Comm_(Comm()),
  Time_(Comm()),
  IsInitialized_(rhs.IsInitialized()),
  IsComputed_(rhs.IsComputed()),
  Condest_(rhs.Condest()),
  LevelOfFill_(rhs.LevelOfFill()),
  H_(0),
  NumInitialize_(rhs.NumInitialize()),
  NumCompute_(rhs.NumCompute()),
  NumApplyInverse_(rhs.NumApplyInverse()),
  InitializeTime_(rhs.InitializeTime()),
  ComputeTime_(rhs.ComputeTime()),
  ApplyInverseTime_(rhs.ApplyInverseTime()),
  Athresh_(rhs.AbsoluteThreshold()),
  Rthresh_(rhs.RelativeThreshold())
{
  H_ = new Epetra_CrsMatrix(rhs.H());
}

//==============================================================================
Ifpack_vIct::~Ifpack_vIct()
{

  if (H_)
    delete H_;

  IsInitialized_ = false;
  IsComputed_ = false;
}

#ifdef HAVE_IFPACK_TEUCHOS
//==========================================================================
int Ifpack_vIct::SetParameters(Teuchos::ParameterList& List)
{

  LevelOfFill_ = List.get("fact: level-of-fill",LevelOfFill());
  Athresh_ = List.get("fact: absolute threshold", Athresh_);
  Rthresh_ = List.get("fact: relative threshold", Rthresh_);

  // set label
  sprintf(Label_, "vICT (fill=%d, athr=%f, rthr=%f)",
	  LevelOfFill(), AbsoluteThreshold(), 
	  RelativeThreshold());
  return(0);
}
#endif

//==========================================================================
int Ifpack_vIct::Initialize()
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
  if (H_)
    delete H_;

  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_.ElapsedTime();

  return(0);
}

//==========================================================================
int Ifpack_vIct::Compute() {

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
    IFPACK_CHK_ERR(-1);

  // insert element Matrix()(0,0)
  Values[0][0] = sqrt(Values[0][0]);
  EPETRA_CHK_ERR(H_->InsertGlobalValues(0,Nnz[0],&(Values[0][0]),
					&(Indices[0][0])));

  if (LevelOfFill() == 0) {

    // ============ //
    // zero fill-in //
    // ============ //

    // now perform the factorization for row `row_i'
    for (int row_i = 1 ; row_i < NumMyRows_ ; ++row_i) {

      // pick col `col_j' (not diagonal)
      for (int j = 0 ; j < Nnz[row_i] - 1 ; ++j) {
	int col_j = Indices[row_i][j];
	double& h_ij = Values[row_i][j];

	// get pointers to row `col_j'
	int* ColIndices;
	double* ColValues;
	int ColNnz;
	H_->ExtractGlobalRowView(col_j, ColNnz, ColValues, ColIndices);

	// look for cross product between row `row_i' and `col_j'
	for (int k = 0 ; k < Nnz[col_j] - 1 ; ++k) {
	  int col_k = ColIndices[k];

	  // col_k is contained in row col_j, look for it
	  // int row `row_i'
	  vector<int>::iterator
	    where = find((Indices[row_i]).begin(), (Indices[row_i]).end(), col_k);
	  if (where != Indices[row_i].end()) {
	    // a common non-zero element exists, add product
	    h_ij -= Values[row_i][*where] * ColValues[k];
	  }
	}
	h_ij /= ColValues[Diag[col_j]];
	Values[row_i][Diag[row_i]] -= h_ij * h_ij;
      }
 
      // adjust diagonal value
      Values[row_i][Diag[row_i]] = sqrt(Values[row_i][Diag[row_i]]);

      // insert row of H
      EPETRA_CHK_ERR(H_->InsertGlobalValues(row_i,Nnz[row_i],
					    &(Values[row_i][0]),
					    &(Indices[row_i][0])));

      Values[row_i].resize(0);
      Indices[row_i].resize(0);
    }

  }
  else {

    // ================ //
    // variable fill-in //
    // ================ //

    vector<double> tmp(NumMyRows_);

    // now perform the factorization for row `row_i'
    for (int row_i = 1 ; row_i < NumMyRows_ ; ++row_i) {

      // number of nonzeros in this row are defined as the nonzeros
      // of the matrix, plus the level of fill (-1, because Nnz
      // contains the diagonal as well, which is always stored).
      int LOF = LevelOfFill() + Nnz[row_i] - 1;

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
	}

	if (h_ij == 0.0 || IFPACK_ABS(h_ij) < AbsoluteThreshold()
	    || IFPACK_ABS(h_ij < DiagonalValue) < RelativeThreshold()) 
	  continue;

	h_ij /= ColValues[ColNnz - 1];
	Element.SetIndex(col_j);
	Element.SetValue(h_ij);
	RowEntries.push_back(Element);
      }

      if (RowEntries.size() > LOF) {
	// sort in ascending order the entries for this row
	sort(RowEntries.begin(),RowEntries.end());
      }
      
      // look for the components that are in the level-of-fill range
      // NOTE: here level-of-fill refers to the number of
      // *off-diagonal* elements. If zero, no off-diagonal elements
      // will be mantained.
      for (int k = 0 ; k < EPETRA_MIN(LOF,RowEntries.size()) ; ++k) {
	int col_k = RowEntries[k].Index();
	if (col_k >= row_i)
	  continue;
	double h_ij = RowEntries[k].Value();
	// skip zero elements
	if (h_ij == 0.0)
	  continue;
	DiagonalValue -= h_ij * h_ij;
	H_->InsertGlobalValues(row_i,1, &h_ij, &col_k);
      }
      // FIXME: use relax here??

      // diagonal element is always inserted
      DiagonalValue = sqrt(DiagonalValue);
      H_->InsertGlobalValues(row_i,1, &DiagonalValue, &row_i);

      Values[row_i].resize(0);
      Indices[row_i].resize(0);
    }
  }

  IFPACK_CHK_ERR(H_->FillComplete());

  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += Time_.ElapsedTime();

  return(0);

}

//=============================================================================
int Ifpack_vIct::ApplyInverse(const Epetra_MultiVector& X, 
			     Epetra_MultiVector& Y) const
{

  Time_.ResetStartTime();

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK_CHK_ERR(-1); // Return error: X and Y not the same size

  if (!IsComputed())
    IFPACK_CHK_ERR(-1); // compute preconditioner first

  // NOTE: I do suppose that X and Y are two different vectors
  EPETRA_CHK_ERR(H_->Solve(false,false,false,X,Y));
  EPETRA_CHK_ERR(H_->Solve(false,true,false,Y,Y));

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_.ElapsedTime();

}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_vIct::Apply(const Epetra_MultiVector& X, 
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
double Ifpack_vIct::Condest(const Ifpack_CondestType CT, 
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
Ifpack_vIct::Print(std::ostream& os) const
{
  os << endl << "*** Ifpack_vIct : " << Label() << endl << endl;
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
    os << "Number of rows of H             = " << H_->NumMyRows() << endl;
    os << "Number of nonzeros of H         = " << H_->NumMyNonzeros() << endl;
    os << "nonzeros / rows                 = " 
       << 1.0 * H_->NumMyNonzeros() / H_->NumMyRows() << endl;
    Epetra_Vector Diagonal(H_->RowMatrixRowMap());
    H_->ExtractDiagonalCopy(Diagonal);
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
