//@HEADER
/*
************************************************************************

		    Komplex: Complex Linear Solver Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#include "Komplex_Operator.hpp"

//==========================================================================
Komplex_Operator::Komplex_Operator(Epetra_DataAccess CV, 
					     Epetra_Operator* Operator,
                                   Komplex_KForms KForm)
 : 
{
}

//==========================================================================
Komplex_Operator::Komplex_Operator(Epetra_DataAccess CV, 
					     Epetra_Operator* Real,
                                   Epetra_Operator* Imag,
					     Komplex_KForms KForm)
 : 
{
}

//==========================================================================
Komplex_Operator::~Komplex_Operator() {
}

//==========================================================================
int Komplex_Operator::SetUseTranspose(bool UseTranspose) {
}

//==========================================================================
int Komplex_Operator::Apply(const Epetra_MultiVector& X, 
				    Epetra_MultiVector& Y) const {
}

//==========================================================================
int Komplex_Operator::ApplyInverse(const Epetra_MultiVector& X, 
					     Epetra_MultiVector& Y) const {
}

//==========================================================================
double Komplex_Operator::NormInf() const {
}

//==========================================================================
const char * Komplex_Operator::Label() const {
}

//==========================================================================
bool Komplex_Operator::UseTranspose() const {
}

//==========================================================================
bool Komplex_Operator::HasNormInf() const {
}

//==========================================================================
const Epetra_Comm & Komplex_Operator::Comm() const {
}

//==========================================================================
const Epetra_Map & Komplex_Operator::OperatorDomainMap() const {
}

//==========================================================================
const Epetra_Map & Komplex_Operator::OperatorRangeMap() const {
}