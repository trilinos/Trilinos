
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER

#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_Util.h"

//=============================================================================
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix()
	: Epetra_Object("Epetra::IntSerialDenseMatrix"),
		CV_(Copy),
		A_Copied_(false),
		M_(0),
		N_(0),
		LDA_(0),
		A_(0)
{
}

//=============================================================================
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix(int NumRows, int NumCols)
  : Epetra_Object("Epetra::IntSerialDenseMatrix"),
		CV_(Copy),
    A_Copied_(false),
    M_(0),
    N_(0),
    LDA_(0),
    A_(0)
{
	if(NumRows < 0)
		throw ReportError("NumRows = " + toString(NumRows) + ". Should be >= 0", -1);
	if(NumCols < 0)
		throw ReportError("NumCols = " + toString(NumCols) + ". Should be >= 0", -1);

	int errorcode = Shape(NumRows, NumCols);
	if(errorcode != 0)
		throw ReportError("Shape returned non-zero (" + toString(errorcode) + ").", -2);
}
//=============================================================================
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix(Epetra_DataAccess CV_in, int* A_in, int lda, 
																												 int NumRows, int NumCols)
  : Epetra_Object("Epetra::IntSerialDenseMatrix"),
		CV_(CV_in),
		A_Copied_(false),
		M_(NumRows),
		N_(NumCols),
		LDA_(lda),
		A_(A_in)
{
	if(A_in == 0)
		throw ReportError("Null pointer passed as A_in parameter.", -3);
	if(NumRows < 0)
		throw ReportError("NumRows = " + toString(NumRows) + ". Should be >= 0", -1);
	if(NumCols < 0)
		throw ReportError("NumCols = " + toString(NumCols) + ". Should be >= 0", -1);
	if(lda < 0)
		throw ReportError("LDA = " + toString(lda) + ". Should be >= 0", -1);

  if(CV_in == Copy) {
    LDA_ = M_;
		const int newsize = LDA_ * N_;
		if(newsize > 0) {
			A_ = new int[newsize];
			CopyMat(A_in, lda, M_, N_, A_, LDA_);
			A_Copied_ = true;
		}
		else {
			A_ = 0;
		}
  }
}
//=============================================================================
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix(const Epetra_IntSerialDenseMatrix& Source)
  : Epetra_Object(Source),
		CV_(Source.CV_),
    A_Copied_(false),
    M_(Source.M_),
    N_(Source.N_),
    LDA_(Source.LDA_),
    A_(Source.A_)
{
	if(CV_ == Copy) {
		LDA_ = M_;
		const int newsize = LDA_ * N_;
		if(newsize > 0) {
			A_ = new int[newsize];
			CopyMat(Source.A_, Source.LDA_, M_, N_, A_, LDA_);
			A_Copied_ = true;
		}
		else {
			A_ = 0;
			A_Copied_ = false;
		}
	}
}
//=============================================================================
int Epetra_IntSerialDenseMatrix::Reshape(int NumRows, int NumCols) {
	if(NumRows < 0 || NumCols < 0)
		return(-1);

	int* A_tmp = 0;
	const int newsize = NumRows * NumCols;

	if(newsize > 0) {
		// Allocate space for new matrix
		A_tmp = new int[newsize];
		for(int k = 0; k < newsize; k++) 
			A_tmp[k] = 0; // Zero out values
		int M_tmp = EPETRA_MIN(M_, NumRows);
		int N_tmp = EPETRA_MIN(N_, NumCols);
		if(A_ != 0) 
			CopyMat(A_, LDA_, M_tmp, N_tmp, A_tmp, NumRows); // Copy principal submatrix of A to new A
	}
  
  CleanupData(); // Get rid of anything that might be already allocated  
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
  A_ = A_tmp; // Set pointer to new A
  A_Copied_ = (newsize>0);
  return(0);
}
//=============================================================================
int Epetra_IntSerialDenseMatrix::Shape(int NumRows, int NumCols) {
	if(NumRows < 0 || NumCols < 0)
		return(-1);

  CleanupData(); // Get rid of anything that might be already allocated
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
	const int newsize = LDA_ * N_;
	if(newsize > 0) {
		A_ = new int[newsize];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for
#endif
		for(int k = 0; k < newsize; k++)
			A_[k] = 0; // Zero out values
		A_Copied_ = true;
	}

  return(0);
}
//=============================================================================
Epetra_IntSerialDenseMatrix::~Epetra_IntSerialDenseMatrix()
{
  CleanupData();
}
//=============================================================================
void Epetra_IntSerialDenseMatrix::CleanupData()
{
  if(A_Copied_)
		delete[] A_; 
	A_ = 0; 
	A_Copied_ = false;
	M_ = 0;
	N_ = 0;
	LDA_ = 0;
}
//=============================================================================
Epetra_IntSerialDenseMatrix& Epetra_IntSerialDenseMatrix::operator = (const Epetra_IntSerialDenseMatrix& Source) {
  if(this == &Source)
		return(*this); // Special case of source same as target
	if((CV_ == View) && (Source.CV_ == View) && (A_ == Source.A_))
		return(*this); // Special case of both are views to same data.

	if(std::strcmp(Label(), Source.Label()) != 0)
		throw ReportError("operator= type mismatch (lhs = " + string(Label()) + 
											", rhs = " + string(Source.Label()) + ").", -5);
	
	if(Source.CV_ == View) {
		if(CV_ == Copy) { // C->V only
			CleanupData();
			CV_ = View;
		}
		M_ = Source.M_; // C->V and V->V
		N_ = Source.N_;
		LDA_ = Source.LDA_;
		A_ = Source.A_;
	}
	else {
		if(CV_ == View) { // V->C
			CV_ = Copy;
			M_ = Source.M_;
			N_ = Source.N_;
			LDA_ = Source.M_;
			const int newsize = LDA_ * N_;
			if(newsize > 0) {
				A_ = new int[newsize];
				A_Copied_ = true;
			}
			else {
				A_ = 0;
				A_Copied_ = false;
			}
		}
		else { // C->C
			if((Source.M_ <= LDA_) && (Source.N_ == N_)) { // we don't need to reallocate
				M_ = Source.M_;
				N_ = Source.N_;
			}
			else { // we need to allocate more space (or less space)
				CleanupData();
				M_ = Source.M_;
				N_ = Source.N_;
				LDA_ = Source.M_;
				const int newsize = LDA_ * N_;
				if(newsize > 0) {
					A_ = new int[newsize];
					A_Copied_ = true;
				}
			}
		}
		CopyMat(Source.A_, Source.LDA_, M_, N_, A_, LDA_); // V->C and C->C
	}
	
  return(*this);
}

//=============================================================================
bool Epetra_IntSerialDenseMatrix::operator==(const Epetra_IntSerialDenseMatrix& rhs) const
{
  if (M_ != rhs.M_ || N_ != rhs.N_) return(false);

  const int* A_tmp = A_;
  const int* rhsA = rhs.A_;

  for(int j=0; j<N_; ++j) {
    int offset = j*LDA_;
    int rhsOffset = j*rhs.LDA_;
    for(int i=0; i<M_; ++i) {
      if (A_tmp[offset+i] != rhsA[rhsOffset+i]) {
	return(false);
      }
    }
  }

  return(true);
}

//=============================================================================
int Epetra_IntSerialDenseMatrix::MakeViewOf(const Epetra_IntSerialDenseMatrix& Source) {
	if(std::strcmp(Label(), Source.Label()) != 0)
		return(-1);

	if(CV_ == Copy) {
		CleanupData();
		CV_ = View;
	}
	M_ = Source.M_;
	N_ = Source.N_;
	LDA_ = Source.LDA_;
	A_ = Source.A_;

	return(0);
}

//=============================================================================
void Epetra_IntSerialDenseMatrix::CopyMat(int* Source, int Source_LDA, int NumRows, int NumCols, 
																					int* Target, int Target_LDA) 
{
  int i, j;
  int* targetPtr = Target;
  int* sourcePtr = 0;
  for(j = 0; j < NumCols; j++) { // for each column
    targetPtr = Target + j * Target_LDA; // set pointers to correct stride
		sourcePtr = Source + j * Source_LDA;
    for(i = 0; i < NumRows; i++) // for each row
			*targetPtr++ = *sourcePtr++; // copy element (i,j) and increment pointer to (i,j+1)
  }
  return;
}

//=============================================================================
int Epetra_IntSerialDenseMatrix::OneNorm() {
	int anorm = 0;
	int* ptr = 0;
	for(int j = 0; j < N_; j++) {
		int sum = 0;
		ptr = A_ + j*LDA_;
		for(int i = 0; i < M_; i++) 
			sum += std::abs(*ptr++);
		anorm = EPETRA_MAX(anorm, sum);
	}
	return(anorm);
}

//=============================================================================
int Epetra_IntSerialDenseMatrix::InfNorm() {	
	int anorm = 0;
	int* ptr = 0;
	// Loop across columns in inner loop.  Most expensive memory access, but 
	// requires no extra storage.
	for(int i = 0; i < M_; i++) {
		int sum = 0;
		ptr = A_ + i;
		for(int j = 0; j < N_; j++) {
			sum += std::abs(*ptr);
			ptr += LDA_;
		}
		anorm = EPETRA_MAX(anorm, sum);
	}
	return(anorm);
}

//=========================================================================
void Epetra_IntSerialDenseMatrix::Print(ostream& os) const {
	if(CV_ == Copy)
		os << "Data access mode: Copy" << endl;
	else
		os << "Data access mode: View" << endl;
	if(A_Copied_)
		os << "A_Copied: yes" << endl;
	else
		os << "A_Copied: no" << endl;
	os << "Rows(M): " << M_ << endl;
	os << "Columns(N): " << N_ << endl;
	os << "LDA: " << LDA_ << endl;
	if(M_ == 0 || N_ == 0)
		os << "(matrix is empty, no values to display)" << endl;
	else
		for(int i = 0; i < M_; i++) {
			for(int j = 0; j < N_; j++){
				os << (*this)(i,j) << " ";
			}
			os << endl;
		}
}

//=========================================================================
int Epetra_IntSerialDenseMatrix::Random() {

	Epetra_Util util;

	for(int j = 0; j < N_; j++) {
		int* arrayPtr = A_ + (j * LDA_);
		for(int i = 0; i < M_; i++) {
			*arrayPtr++ = util.RandomInt();
		}
	}
	
	return(0);
}
