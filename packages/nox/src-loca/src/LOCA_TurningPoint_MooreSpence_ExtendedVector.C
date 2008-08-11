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

#include "LOCA_TurningPoint_MooreSpence_ExtendedVector.H"  // Class definition
#include "LOCA_TurningPoint_MooreSpence_ExtendedMultiVector.H"

LOCA::TurningPoint::MooreSpence::ExtendedVector::ExtendedVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& xVec,
		    const NOX::Abstract::Vector& nullVec,
		    double bifParam) :
  LOCA::Extended::Vector(global_data,2,1)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, bifParam);
}

LOCA::TurningPoint::MooreSpence::ExtendedVector::ExtendedVector(
                const LOCA::TurningPoint::MooreSpence::ExtendedVector& source,
		NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::TurningPoint::MooreSpence::ExtendedVector::~ExtendedVector()
{
}

NOX::Abstract::Vector& 
LOCA::TurningPoint::MooreSpence::ExtendedVector::operator=(
					      const NOX::Abstract::Vector& y)
{
  operator=(dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(y));
  return *this;
}

LOCA::Extended::Vector& 
LOCA::TurningPoint::MooreSpence::ExtendedVector::operator=(
					     const LOCA::Extended::Vector& y)
{
  operator=(dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(y));
  return *this;
}

LOCA::TurningPoint::MooreSpence::ExtendedVector& 
LOCA::TurningPoint::MooreSpence::ExtendedVector::operator=(
                     const LOCA::TurningPoint::MooreSpence::ExtendedVector& y)
{ 
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::clone(
						    NOX::CopyType type) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedVector(*this, 
								     type));
}

void 
LOCA::TurningPoint::MooreSpence::ExtendedVector::setVec(
					const NOX::Abstract::Vector& xVec,
					const NOX::Abstract::Vector& nullVec,
					double bifPar)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, bifPar);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::getXVec() const
{
  return getVector(0);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::getNullVec() const
{
  return getVector(1);
}

double 
LOCA::TurningPoint::MooreSpence::ExtendedVector::getBifParam() const
{
  return getScalar(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::getXVec()
{
  return getVector(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::getNullVec()
{
  return getVector(1);
}

double& 
LOCA::TurningPoint::MooreSpence::ExtendedVector::getBifParam()
{
  return getScalar(0);
}

LOCA::TurningPoint::MooreSpence::ExtendedVector::ExtendedVector(
		  const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  LOCA::Extended::Vector(global_data,2,1)
{
}

Teuchos::RCP<LOCA::Extended::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::generateMultiVector(
							int nColumns, 
							int nVectorRows, 
							int nScalarRows) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedMultiVector(
								    globalData,
								    nColumns));
}
