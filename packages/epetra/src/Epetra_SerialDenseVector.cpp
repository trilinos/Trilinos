
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
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
/*double& Epetra_SerialDenseVector::operator() (int Index)  {
  if (Index >= M_ || Index < 0)
		throw ReportError("Index = " +toString(Index) + " Out of Range 0 - " + toString(M_-1), -1); 
  return(A_[Index]); 
} 
//=========================================================================
const double& Epetra_SerialDenseVector::operator() (int Index) const  { 
  if (Index >= M_ || Index < 0)  
		throw ReportError("Index = " +toString(Index) + " Out of Range 0 - " + toString(M_-1), -1); 
   return(A_[Index]); 
}
//=========================================================================
const double& Epetra_SerialDenseVector::operator [] (int Index) const  { 
   return(A_[Index]); 
} 
//=========================================================================
double& Epetra_SerialDenseVector::operator [] (int Index)  { 
  if (Index >= M_ || Index < 0)  
		throw ReportError("Index = " +toString(Index) + " Out of Range 0 - " + toString(M_-1), -1); 
   return(A_[Index]); 
}*/

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
