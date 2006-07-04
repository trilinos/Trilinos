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

#include "LOCA_Abstract_Group.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Abstract::Group::Group(
		     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data)
  : LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup(Teuchos::rcp(new LOCA::DerivUtils(global_data))),
    globalData(global_data)
{
}

LOCA::Abstract::Group::Group(
		     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		     const Teuchos::RefCountPtr<LOCA::DerivUtils>& d)
  : LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup(d),
    globalData(global_data)
{
}

// LOCA::Abstract::Group::Group(Teuchos::ParameterList& params, 
// 			     const LOCA::DerivUtils& d)
//   : LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup(d)
// {
// }

LOCA::Abstract::Group::Group(const LOCA::Abstract::Group& source, 
			     NOX::CopyType type) : 
  LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup(source, type),
    globalData(source.globalData)
{
}


LOCA::Abstract::Group::~Group() 
{
}

// NOX::Abstract::Group::ReturnType 
// LOCA::Abstract::Group::augmentJacobianForHomotopy(double conParamValue)
// {
//   return NOX::Abstract::Group::NotDefined;
// }

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyShiftedMatrixInverse(
					   Teuchos::ParameterList& params, 
					   const NOX::Abstract::Vector& input,
					   NOX::Abstract::Vector& result,
					   double shift)
{
  globalData->locaErrorCheck->throwError(
			   "LOCA::Abstract::Group::applyShiftedMatrixInverse",
			   "Not implemented for group");
  return NOX::Abstract::Group::NotDefined;
}

// NOX::Abstract::Group::ReturnType
// LOCA::Abstract::Group::applyComplexInverse(
// 			       Teuchos::ParameterList& params,
// 			       const NOX::Abstract::Vector& input_real,
// 			       const NOX::Abstract::Vector& input_imag,
// 			       double frequency,
// 			       NOX::Abstract::Vector& result_real,
// 			       NOX::Abstract::Vector& result_imag) const
// {
//   globalData->locaErrorCheck->throwError("LOCA::Abstract::Group::applyComplexInverse",
// 			       "No mass matrix defined for group");
//   return NOX::Abstract::Group::NotDefined;
// }

void
LOCA::Abstract::Group::copy(const NOX::Abstract::Group& src)
{
  const LOCA::Abstract::Group& source = 
    dynamic_cast<const LOCA::Abstract::Group&>(src);

  // Copy parent classes
  LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::copy(src);
  globalData = source.globalData;
}

void
LOCA::Abstract::Group::setParamsMulti(const vector<int>& paramIDs, 
		const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  for (unsigned int i=0; i<paramIDs.size(); i++)
    setParam(paramIDs[i], vals(i,0));
}

void
LOCA::Abstract::Group::notifyCompletedStep()
{
}

NOX::Abstract::Group&
LOCA::Abstract::Group::operator=(const NOX::Abstract::Group& source)
{

  copy(source);
  return *this;
}
