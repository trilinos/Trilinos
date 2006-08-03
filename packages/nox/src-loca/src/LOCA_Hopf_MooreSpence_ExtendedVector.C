// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_Hopf_MooreSpence_ExtendedVector.H"  // Class definition
#include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H"

LOCA::Hopf::MooreSpence::ExtendedVector::ExtendedVector(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& xVec,
		    const NOX::Abstract::Vector& realEigenVec,
		    const NOX::Abstract::Vector& imagEigenVec,
		    double frequency,
		    double bifParam) :
  LOCA::Extended::Vector(global_data,3,2)
{
  setVector(0, xVec);
  setVector(1, realEigenVec);
  setVector(2, imagEigenVec);
  setScalar(0, frequency);
  setScalar(1, bifParam);
}

LOCA::Hopf::MooreSpence::ExtendedVector::ExtendedVector(
                    const LOCA::Hopf::MooreSpence::ExtendedVector& source,
		    NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::Hopf::MooreSpence::ExtendedVector::~ExtendedVector()
{
}

NOX::Abstract::Vector& 
LOCA::Hopf::MooreSpence::ExtendedVector::operator=(
					      const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedVector&>(y));
}

LOCA::Extended::Vector& 
LOCA::Hopf::MooreSpence::ExtendedVector::operator=(
					      const LOCA::Extended::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedVector&>(y));
}

LOCA::Hopf::MooreSpence::ExtendedVector& 
LOCA::Hopf::MooreSpence::ExtendedVector::operator=(
                         const LOCA::Hopf::MooreSpence::ExtendedVector& y)
{ 
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

Teuchos::RefCountPtr<NOX::Abstract::Vector> 
LOCA::Hopf::MooreSpence::ExtendedVector::clone(NOX::CopyType type) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedVector(*this, type));
}

void 
LOCA::Hopf::MooreSpence::ExtendedVector::setVec(
			       const NOX::Abstract::Vector& xVec,
			       const NOX::Abstract::Vector& realEigenVec,
			       const NOX::Abstract::Vector& imagEigenVec,
			       double frequency,
			       double bifPar)
{
  setVector(0, xVec);
  setVector(1, realEigenVec);
  setVector(2, imagEigenVec);
  setScalar(0, frequency);
  setScalar(1, bifPar);
}

Teuchos::RefCountPtr<const NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getXVec() const
{
  return getVector(0);
}

Teuchos::RefCountPtr<const NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getRealEigenVec() const
{
  return getVector(1);
}

Teuchos::RefCountPtr<const NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getImagEigenVec() const
{
  return getVector(2);
}

double 
LOCA::Hopf::MooreSpence::ExtendedVector::getFrequency() const
{
  return getScalar(0);
}

double 
LOCA::Hopf::MooreSpence::ExtendedVector::getBifParam() const
{
  return getScalar(1);
}

Teuchos::RefCountPtr<NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getXVec()
{
  return getVector(0);
}

Teuchos::RefCountPtr<NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getRealEigenVec()
{
  return getVector(1);
}

Teuchos::RefCountPtr<NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getImagEigenVec()
{
  return getVector(2);
}

double& 
LOCA::Hopf::MooreSpence::ExtendedVector::getFrequency()
{
  return getScalar(0);
}

double& 
LOCA::Hopf::MooreSpence::ExtendedVector::getBifParam()
{
  return getScalar(1);
}

LOCA::Hopf::MooreSpence::ExtendedVector::ExtendedVector(
		  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) :
  LOCA::Extended::Vector(global_data,3,2)
{
}

Teuchos::RefCountPtr<LOCA::Extended::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedVector::generateMultiVector(
							int nColumns, 
							int nVectorRows, 
							int nScalarRows) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedMultiVector(globalData,
								  nColumns));
}
