
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

#include "Epetra_SerialDenseVector.h"
//=============================================================================
Epetra_SerialDenseVector::Epetra_SerialDenseVector()
  : Epetra_SerialDenseMatrix()
{
	SetLabel("Epetra::SerialDenseVector");
}

//=============================================================================
Epetra_SerialDenseVector::Epetra_SerialDenseVector(int length)
  : Epetra_SerialDenseMatrix(length, 1)
{
	SetLabel("Epetra::SerialDenseVector");
}

//=============================================================================
Epetra_SerialDenseVector::Epetra_SerialDenseVector(Epetra_DataAccess CV_in, double *values, int length)
  : Epetra_SerialDenseMatrix(CV_in, values, length, length, 1)
{
	SetLabel("Epetra::SerialDenseVector");	
}

//=============================================================================
Epetra_SerialDenseVector::Epetra_SerialDenseVector(const Epetra_SerialDenseVector& Source)
  : Epetra_SerialDenseMatrix(Source)
{}

//=============================================================================
Epetra_SerialDenseVector::~Epetra_SerialDenseVector()
{}

//=========================================================================
Epetra_SerialDenseVector& Epetra_SerialDenseVector::operator = (const Epetra_SerialDenseVector& Source) {
	Epetra_SerialDenseMatrix::operator=(Source); // call this->Epetra_SerialDenseMatrix::operator =
	return(*this);
}

//=========================================================================
int Epetra_SerialDenseVector::Random() {
	int errorcode = Epetra_SerialDenseMatrix::Random();
	return(errorcode);
}

//=========================================================================
double  Epetra_SerialDenseVector::Dot(const Epetra_SerialDenseVector & x) const {
  
#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
  if (Length()!=x.Length()) 
    throw ReportError("Length of this object = " + 
		      toString(Length()) + " is not equal to length of x = "  + toString(x.Length()), -1);
#endif

  // dot-product of this and x.
    
  double result = DOT(Length(), Values(), x.Values());
  
  UpdateFlops(2*Length());

  return(result);
}

//=========================================================================
double  Epetra_SerialDenseVector::Norm1() const {
  
  // 1-norm of vector
    
  double result = ASUM(Length(), Values());
  
  UpdateFlops(2*Length());

  return(result);
}

//=========================================================================
double  Epetra_SerialDenseVector::Norm2() const {
  
  // 2-norm of vector
    
  double result = NRM2(Length(), Values());
  
  UpdateFlops(2*Length());

  return(result);
}
//=========================================================================
double  Epetra_SerialDenseVector::NormInf() const {
  
  // Inf-norm of vector
  double result = 0.0;
  int j = IAMAX(Length(), Values()); // Location of max (-1) if length zero

  if (j>-1) result = std::abs( (*this)[j]);
  
  // UpdateFlops(2*Length()); // Technically there are no FLOPS

  return(result);
}
//=========================================================================
void Epetra_SerialDenseVector::Print(ostream& os) const {
	if(CV_ == Copy)
		os << "Data access mode: Copy" << endl;
	else
		os << "Data access mode: View" << endl;
	if(A_Copied_)
		os << "A_Copied: yes" << endl;
	else
		os << "A_Copied: no" << endl;
	os << "Length(M): " << M_ << endl;
	if(M_ == 0)
		os << "(vector is empty, no values to display)";
	else
		for(int i = 0; i < M_; i++)
      os << (*this)(i) << " ";
	os << endl;
}
