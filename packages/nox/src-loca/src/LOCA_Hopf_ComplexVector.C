// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_Hopf_ComplexVector.H"  // Class definition
#include "LOCA_Hopf_ComplexMultiVector.H"

LOCA::Hopf::ComplexVector::ComplexVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& realVec,
		    const NOX::Abstract::Vector& imagVec) :
  LOCA::Extended::Vector(global_data,2,0)
{
  setVector(0, realVec);
  setVector(1, imagVec);
}

LOCA::Hopf::ComplexVector::ComplexVector(
				     const LOCA::Hopf::ComplexVector& source,
				     NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::Hopf::ComplexVector::~ComplexVector()
{
}

NOX::Abstract::Vector& 
LOCA::Hopf::ComplexVector::operator=(const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Hopf::ComplexVector&>(y));
}

LOCA::Extended::Vector& 
LOCA::Hopf::ComplexVector::operator=(const LOCA::Extended::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Hopf::ComplexVector&>(y));
}

LOCA::Hopf::ComplexVector& 
LOCA::Hopf::ComplexVector::operator=(const LOCA::Hopf::ComplexVector& y)
{ 
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector> 
LOCA::Hopf::ComplexVector::clone(NOX::CopyType type) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::ComplexVector(*this, type));
}

void 
LOCA::Hopf::ComplexVector::setVec(const NOX::Abstract::Vector& realVec,
				  const NOX::Abstract::Vector& imagVec)
{
  setVector(0, realVec);
  setVector(1, imagVec);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::ComplexVector::getRealVec() const
{
  return getVector(0);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::ComplexVector::getImagVec() const
{
  return getVector(1);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Hopf::ComplexVector::getRealVec()
{
  return getVector(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Hopf::ComplexVector::getImagVec()
{
  return getVector(1);
}

LOCA::Hopf::ComplexVector::ComplexVector(
		  const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  LOCA::Extended::Vector(global_data,2,0)
{
}

Teuchos::RCP<LOCA::Extended::MultiVector>
LOCA::Hopf::ComplexVector::generateMultiVector(int nColumns, 
					       int nVectorRows, 
					       int nScalarRows) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::ComplexMultiVector(globalData, nColumns));
}
