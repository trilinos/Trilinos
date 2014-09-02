
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

#include "Ifpack_SerialTriDiMatrix.h"
#include "Epetra_Util.h"

//=============================================================================
Ifpack_SerialTriDiMatrix::Ifpack_SerialTriDiMatrix(bool set_object_label)
  : Epetra_CompObject(),
    Epetra_Object(-1, false),
    N_(0),
    A_Copied_(false),
    CV_(Copy),
    A_(0),
    UseTranspose_(false)
{
  if (set_object_label) {
    SetLabel("Epetra::SerialTriDiMatrix");
  }
}

//=============================================================================
Ifpack_SerialTriDiMatrix::Ifpack_SerialTriDiMatrix(int NumRowCol,
                                                   bool set_object_label)
  : Epetra_CompObject(),
    Epetra_Object(-1, false),
    N_(NumRowCol),
    A_Copied_(false),
    CV_(Copy),
    A_(0),
    UseTranspose_(false)
{
  if (set_object_label) {
    SetLabel("Epetra::SerialTriDiMatrix");
  }
  if(NumRowCol < 0)
	throw ReportError("NumRows = " + toString(NumRowCol) + ". Should be >= 0", -1);

  int errorcode = Shape(NumRowCol);
  if(errorcode != 0)
    throw ReportError("Shape returned non-zero value", errorcode);
}

//=============================================================================
Ifpack_SerialTriDiMatrix::Ifpack_SerialTriDiMatrix(Epetra_DataAccess CV_in, double* A_in,
                                                   int NumRowCol,
                                                   bool set_object_label)
  : Epetra_CompObject(),
    Epetra_Object(-1, false),
    N_(NumRowCol),
    A_Copied_(false),
    CV_(CV_in),
    A_(A_in),
    UseTranspose_(false)
{
  if (set_object_label) {
    SetLabel("Epetra::SerialTriDiMatrix");
  }
  if(A_in == 0)
	throw ReportError("Null pointer passed as A parameter.", -3);
  if(NumRowCol < 0)
    throw ReportError("NumRowCol = " + toString(NumRowCol) + ". Should be >= 0", -1);

  if (CV_in == Copy) {
    const int newsize = (N_ == 1) ? 1 : 4*(N_-1);
    if (newsize > 0) {
      A_ = new double[newsize];
      CopyMat(A_in, N_, A_, N_);
      A_Copied_ = true;
      DL_ = A_;
      D_  = DL_+(N_-1);
      DU_ = D_ + N_;
      DU2_ = DU_ + (N_-1);
    }
    else {
      A_ = 0;
    }
  }
}
//=============================================================================
Ifpack_SerialTriDiMatrix::Ifpack_SerialTriDiMatrix(const Ifpack_SerialTriDiMatrix& Source)
  : Epetra_CompObject(Source),
    N_(Source.N_),
    LDA_(Source.LDA_),
    A_Copied_(false),
    CV_(Source.CV_),
    UseTranspose_(false)
{
  SetLabel(Source.Label());
  if(CV_ == Copy) {
    const int newsize = 4* (N_-1);
    if(newsize > 0) {
      A_ = new double[newsize];
      CopyMat(Source.A_, Source.N() , A_, N_);
      A_Copied_ = true;
      DL_ = A_;
      D_  = DL_+(N_-1);
      DU_ = D_ + N_;
      DU2_ = DU_ + (N_-1);
    }
    else {
      A_ = 0;
    }
  }
}

//=============================================================================
int Ifpack_SerialTriDiMatrix::Reshape(int NumRows, int NumCols) {
	if(NumRows < 0 || NumCols < 0)
	  return(-1);
	if(NumRows != NumCols)
	  return(-1);

	double* A_tmp = 0;
	const int newsize = (N_ == 1)? 1 : 4*(N_-1);
	if(newsize > 0) {
		// Allocate space for new matrix
		A_tmp = new double[newsize];
		for (int k = 0; k < newsize; k++)
			A_tmp[k] = 0.0; // Zero out values
		if (A_ != 0)
		  CopyMat(A_, N_, A_tmp, NumRows); // Copy principal submatrix of A to new A
		
  }
  CleanupData(); // Get rid of anything that might be already allocated
  N_ = NumCols;
  if(newsize > 0) {
    A_ = A_tmp; // Set pointer to new A
    A_Copied_ = true;
  }

  DL_ = A_;
  D_ = A_+(N_-1);
  DU_ = D_+N_;
  DU2_ = DU_+(N_-1);

  return(0);
}
//=============================================================================
int Ifpack_SerialTriDiMatrix::Shape(int NumRowCol) {
	if(NumRowCol < 0 || NumRowCol < 0)
		return(-1);

  CleanupData(); // Get rid of anything that might be already allocated
  N_ = NumRowCol;
  LDA_ = N_;
  const int newsize = (N_ == 1)? 1 : 4*(N_-1);
  if(newsize > 0) {
    A_ = new double[newsize];
    for (int k = 0; k < newsize; k++)
      A_[k] = 0.0; // Zero out values
    A_Copied_ = true;
  }
  // set the pointers
  DL_ = A_;
  D_ = A_+(N_-1);
  DU_ = D_+N_;
  DU2_ = DU_+(N_-1);

  return(0);
}
//=============================================================================
Ifpack_SerialTriDiMatrix::~Ifpack_SerialTriDiMatrix()
{
  CleanupData();
}
//=============================================================================
void Ifpack_SerialTriDiMatrix::CleanupData()
{
  if (A_Copied_)
    delete [] A_;
	A_ = DL_ = D_ = DU_ = DU2_ =  0;
	A_Copied_ = false;
	N_ = 0;
}
//=============================================================================
Ifpack_SerialTriDiMatrix& Ifpack_SerialTriDiMatrix::operator = (const Ifpack_SerialTriDiMatrix& Source) {
  if(this == &Source)
		return(*this); // Special case of source same as target
	if((CV_ == View) && (Source.CV_ == View) && (A_ == Source.A_))
		return(*this); // Special case of both are views to same data.

	if(std::strcmp(Label(), Source.Label()) != 0)
		throw ReportError("operator= type mismatch (lhs = " + std::string(Label()) +
      ", rhs = " + std::string(Source.Label()) + ").", -5);
	
	if(Source.CV_ == View) {
		if(CV_ == Copy) { // C->V only
			CleanupData();
			CV_ = View;
		}
		N_ = Source.N_;
		A_ = Source.A_;
		DL_ = Source.DL_;
		D_ = Source.D_;
		DU_ = Source.DU_;
		DU2_ = Source.DU2_;
	}
	else {
		if(CV_ == View) { // V->C
			CV_ = Copy;
			N_ = Source.N_;
			const int newsize = 4*N_ - 4;
			if(newsize > 0) {
				A_ = new double[newsize];
				DL_ = A_;
				D_ = A_+(N_-1);
				DU_ = D_+N_;
				DU2_ = DU_+(N_-1);
				A_Copied_ = true;
			}
			else {
				A_ = 0;
				DL_ = D_ = DU_ = DU2_ = 0;
				A_Copied_ = false;
			}
		}
		else { // C->C
			if(Source.N_ == N_) { // we don't need to reallocate
				N_ = Source.N_;
			}
			else { // we need to allocate more space (or less space)
				CleanupData();
				N_ = Source.N_;
				const int newsize = (N_ == 1)? 1 : 4*(N_-1);
				if(newsize > 0) {
					A_ = new double[newsize];
					DL_ = A_;
					D_ = A_+(N_-1);
					DU_ = D_+N_;
					DU2_ = DU_+(N_-1);
					A_Copied_ = true;
				}
			}
		}
		CopyMat(Source.A_, Source.N(), A_, N_); // V->C and C->C
	}
	
  return(*this);
}


//=============================================================================
bool Ifpack_SerialTriDiMatrix::operator==(const Ifpack_SerialTriDiMatrix& rhs) const
{
  if (N_ != rhs.N_) return(false);

  const double* A_tmp = A_;
  const double* rhsA = rhs.A_;

  const int size = (N_ == 1)? 1 : 4*(N_-1);

  for(int j=0; j<size; ++j) {
    if (std::abs(A_tmp[j] - rhsA[j]) > Epetra_MinDouble) {
      return(false);
    }
  }

  return(true);
}

//=============================================================================
Ifpack_SerialTriDiMatrix& Ifpack_SerialTriDiMatrix::operator+= ( const Ifpack_SerialTriDiMatrix & Source) {
    if (N() != Source.N())
		throw ReportError("Column dimension of source = " + toString(Source.N()) +
											" is different than column dimension of target = " + toString(N()), -2);

    CopyMat(Source.A(), Source.N(), A(), N(),true);
  return(*this);
}
//=============================================================================
void Ifpack_SerialTriDiMatrix::CopyMat(const double* Source,
				       int nrowcol,
                                       double* Target,
				       int tN,
                                       bool add)
{
  int lmax = EPETRA_MIN(nrowcol,tN);
  if (add) {
    // do this in 4 passes
    for(int j=0; j<lmax; ++j) {      
      Target[(tN-1)+j] += Source[(nrowcol-1)+j];  //D      
      if(j<tN-1) {
	Target[j] += Source[j];  // DL
	Target[(tN-1)+tN + j] += Source[(nrowcol-1)+ nrowcol + j]; // DU
      }
      if(j<tN-2) Target[(tN-1)*2 + tN + j] += Source[ (nrowcol-1)*2 +nrowcol  + j]; // DU2
    }
  }
  else {  
    for(int j=0; j<lmax; ++j) {
      Target[(tN-1)+j] = Source[(nrowcol-1)+j];  //D      
      if(j<tN-1) {
	Target[j] = Source[j];  // DL
	Target[(tN-1)+tN + j] = Source[(nrowcol-1)+ nrowcol + j]; // DU
      }
      if(j<tN-2) Target[(tN-1)*2 + tN + j] = Source[ (nrowcol-1)*2 +nrowcol  + j]; // DU2
    }
  }
  return;
}
//=============================================================================
double Ifpack_SerialTriDiMatrix::NormOne() const {
  int i;
  double anorm = 0.0;
  double * ptr  = A_;
  double sum=0.0;
  
  const int size = (N_ == 1)? 1 : 4*(N_-1);

  for (i=0; i<size; i++) sum += std::abs(*ptr++);
  
  anorm = EPETRA_MAX(anorm, sum);
  UpdateFlops((double)size );
  return(anorm);
}
//=============================================================================
double Ifpack_SerialTriDiMatrix::NormInf() const {
  return NormOne();

}
//=============================================================================
int Ifpack_SerialTriDiMatrix::Scale(double ScalarA) {

  int i;

  double * ptr = A_;

  const int size = (N_ == 1)? 1 : 4*(N_-1);

  for (i=0; i<size ; i++) { *ptr = ScalarA * (*ptr); ptr++; }
 
  UpdateFlops((double)N_*(double)N_);
  return(0);

}

// //=========================================================================
int Ifpack_SerialTriDiMatrix::Multiply (char TransA, char TransB, double ScalarAB,
					const Ifpack_SerialTriDiMatrix& A_in,
					const Ifpack_SerialTriDiMatrix& B,
					double ScalarThis ) {
  // // Check for compatible dimensions

  // if (TransA!='T' && TransA!='N') EPETRA_CHK_ERR(-2); // Return error
  // if (TransB!='T' && TransB!='N') EPETRA_CHK_ERR(-3);

  // // int A_nrows = (TransA=='T') ? A_in.N() : A_in.M();
  // // int A_ncols = (TransA=='T') ? A_in.M() : A_in.N();
  // // int B_nrows = (TransB=='T') ? B.N() : B.M();
  // // int B_ncols = (TransB=='T') ? B.M() : B.N();

  // // if (M_        != A_nrows     ||
  // //     A_ncols   != B_nrows     ||
  // //     N_        != B_ncols  ) EPETRA_CHK_ERR(-1); // Return error


  // // Call GEMM function
  // GEMM(TransA, TransB, M_, N_, A_ncols, ScalarAB, A_in.A(), A_in.LDA(),
  //      B.A(), B.LDA(), ScalarThis, A_, LDA_);
  // long int nflops = 2*M_;
  // nflops *= N_;
  // nflops *= A_ncols;
  // if (ScalarAB != 1.0) nflops += M_*N_;
  // if (ScalarThis != 0.0) nflops += M_*N_;
  // UpdateFlops((double)nflops);

  // return(0);

  return(-1);

}

//=========================================================================
// int  Ifpack_SerialTriDiMatrix::Multiply (bool transA,
//                                          const Ifpack_SerialTriDiMatrix& x,
//                                          Ifpack_SerialTriDiMatrix& y)
// {
//   int A_nrows = M();
//   int x_nrows = x.M();
//   int y_nrows = y.M();
//   int A_ncols = N();
//   int x_ncols = x.N();
//   int y_ncols = y.N();

//   if (transA) {
//     if (x_nrows != A_nrows) EPETRA_CHK_ERR(-1);
//     if (y_ncols != x_ncols || y_nrows != A_ncols) y.Reshape(A_ncols, x_ncols);
//   }
//   else {
//     if (x_nrows != A_ncols) EPETRA_CHK_ERR(-1);
//     if (y_ncols != x_ncols || y_nrows != A_nrows) y.Reshape(A_nrows, x_ncols);
//   }

//   double scalar0 = 0.0;
//   double scalar1 = 1.0;

//   int err = 0;
//   if (transA) {
//     err = y.Multiply('T', 'N', scalar1, *this, x, scalar0);
//   }
//   else {
//     err = y.Multiply('N', 'N', scalar1, *this, x, scalar0);
//   }
//   // FIXME (mfh 06 Mar 2013) Why aren't we returning err instead of 0?
//   // In any case, I'm going to silence the unused value compiler
//   // warning for now.  I'm not changing the return value because I
//   // don't want to break any downstream code that depends on this
//   // method always returning 0.
//   (void) err;

//   return(0);
// }


//=========================================================================
void Ifpack_SerialTriDiMatrix::Print(std::ostream& os) const {
	os << std::endl;
	if(CV_ == Copy)
		os << "Data access mode: Copy" << std::endl;
	else
		os << "Data access mode: View" << std::endl;
	if(A_Copied_)
		os << "A_Copied: yes" << std::endl;
	else
		os << "A_Copied: no" << std::endl;
	os << "Rows and Columns(N): " << N_ << std::endl;
	if(N_ == 0)
		os << "(matrix is empty, no values to display)" << std::endl;
	else
	  {
	    os << "DL: "<<std::endl;
	    for(int i=0;i<N_-1;++i) 
	      os << DL_[i]<<" ";	    
	    os << std::endl;
	    os << "D: "<<std::endl;
	    for(int i=0;i<N_;++i) 
	      os << D_[i]<<" ";	    
	    os << std::endl;
	    os << "DU: "<<std::endl;
	    for(int i=0;i<N_-1;++i) 
	      os << DU_[i]<<" ";	    
	    os << std::endl;
	    os << "DU2: "<<std::endl;
	    for(int i=0;i<N_-2;++i) 
	      os << DU2_[i]<<" ";	    
	    os << std::endl;
	  }
}
//=========================================================================
int Ifpack_SerialTriDiMatrix::Random() {
  
  Epetra_Util util;
  double* arrayPtr = A_;
  const int size = (N_ == 1)? 1 : 4*(N_-1);
  for(int j = 0; j < size ; j++) {	  
    *arrayPtr++ = util.RandomDouble();	
  }
  
  return(0);
}

// //=========================================================================
// int Ifpack_SerialTriDiMatrix::Apply(const Ifpack_SerialTriDiMatrix& X,
// 				    Ifpack_SerialTriDiMatrix& Y) {
//   return Multiply(UseTranspose(), X, Y);
// }
