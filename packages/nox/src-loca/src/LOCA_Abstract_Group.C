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
#include "NOX_Parameter_List.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Abstract::Group::Group(
		     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data)
  : LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup(Teuchos::rcp(new LOCA::DerivUtils(global_data))),
    globalData(global_data)
{
}

LOCA::Abstract::Group::Group(
		     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		     const Teuchos::RefCountPtr<LOCA::DerivUtils>& d)
  : LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup(d),
    globalData(global_data)
{
}

// LOCA::Abstract::Group::Group(NOX::Parameter::List& params, 
// 			     const LOCA::DerivUtils& d)
//   : LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup(d)
// {
// }

LOCA::Abstract::Group::Group(const LOCA::Abstract::Group& source, 
			     NOX::CopyType type)
  : LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup(source, type),
    globalData(source.globalData)
{
}


LOCA::Abstract::Group::~Group() 
{
}

// LOCA::Continuation::AbstractGroup&
// LOCA::Abstract::Group::operator=(
// 			    const LOCA::Continuation::AbstractGroup& source)
// {
//   return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
// }

// LOCA::Continuation::FiniteDifferenceGroup&
// LOCA::Abstract::Group::operator=(
// 		      const LOCA::Continuation::FiniteDifferenceGroup& source)
// {
//   return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
// }

// LOCA::Bifurcation::TPBord::AbstractGroup&
// LOCA::Abstract::Group::operator=(
// 		     const LOCA::Bifurcation::TPBord::AbstractGroup& source)
// {
//   return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
// }

// LOCA::Bifurcation::TPBord::FiniteDifferenceGroup&
// LOCA::Abstract::Group::operator=(
// 	      const LOCA::Bifurcation::TPBord::FiniteDifferenceGroup& source)
// {
//   return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
// }

// LOCA::Bifurcation::TPBord::SingularSolveGroup&
// LOCA::Abstract::Group::operator=(
// 		   const LOCA::Bifurcation::TPBord::SingularSolveGroup& source)
// {
//   return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
// }

LOCA::TimeDependent::AbstractGroup&
LOCA::Abstract::Group::operator=(
			    const LOCA::TimeDependent::AbstractGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

// LOCA::Bifurcation::HopfBord::AbstractGroup&
// LOCA::Abstract::Group::operator=(
// 	       const LOCA::Bifurcation::HopfBord::AbstractGroup& source)
// {
//   return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
// }

// LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup&
// LOCA::Abstract::Group::operator=(
// 	       const LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup& source)
// {
//   return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
// }

// LOCA::Homotopy::AbstractGroup&
// LOCA::Abstract::Group::operator=(
// 			    const LOCA::Homotopy::AbstractGroup& source)
// {
//   return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
// }

LOCA::Abstract::Group&
LOCA::Abstract::Group::operator=(const LOCA::Abstract::Group& source)
{

  // Copy parent classes
//   LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup::operator=(source);
//   LOCA::Bifurcation::TPBord::SingularSolveGroup::operator=(source);
  LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::operator=(source);
  globalData = source.globalData;
  return *this;
}

// NOX::Abstract::Group::ReturnType 
// LOCA::Abstract::Group::augmentJacobianForHomotopy(double conParamValue)
// {
//   return NOX::Abstract::Group::NotDefined;
// }

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyShiftedMatrixInverse(
					   NOX::Parameter::List& params, 
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
// 			       NOX::Parameter::List& params,
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

// NOX::Abstract::Group::ReturnType
// LOCA::Abstract::Group::applyBorderedJacobianInverse(bool trans,
// 				     NOX::Parameter::List& params,
// 				     const NOX::Abstract::Vector& a,
// 				     const NOX::Abstract::Vector& b,
// 				     const NOX::Abstract::Vector& vInput,
// 				     double sInput,
// 				     NOX::Abstract::Vector& vResult,
// 				     double& sResult) const
// {
//   globalData->locaErrorCheck->throwError(
// 		       "LOCA::Abstract::Group::applyBorderedJacobianInverse",
// 		       "Not defined for group");
//   return NOX::Abstract::Group::NotDefined;
// }

LOCA::MultiContinuation::AbstractGroup&
LOCA::Abstract::Group::operator=(
			 const LOCA::MultiContinuation::AbstractGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

void
LOCA::Abstract::Group::setParamsMulti(const vector<int>& paramIDs, 
		const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  for (unsigned int i=0; i<paramIDs.size(); i++)
    setParam(paramIDs[i], vals(i,0));
}

LOCA::MultiContinuation::FiniteDifferenceGroup&
LOCA::Abstract::Group::operator=(
		 const LOCA::MultiContinuation::FiniteDifferenceGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::TurningPoint::MooreSpence::AbstractGroup&
LOCA::Abstract::Group::operator=(
	 const LOCA::TurningPoint::MooreSpence::AbstractGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup&
LOCA::Abstract::Group::operator=(
	 const LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}
