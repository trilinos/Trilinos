
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
Epetra_SerialDenseVector::Epetra_SerialDenseVector(int Length)
  : Epetra_SerialDenseMatrix(Length, 1)
{
	SetLabel("Epetra::SerialDenseVector");
}

//=============================================================================
Epetra_SerialDenseVector::Epetra_SerialDenseVector(Epetra_DataAccess CV, double *Values, int Length)
  : Epetra_SerialDenseMatrix(CV, Values, Length, Length, 1)
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
