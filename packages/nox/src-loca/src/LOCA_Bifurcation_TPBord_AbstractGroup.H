// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_BIFURCATION_TPBORD_ABSTRACTGROUP_H
#define LOCA_BIFURCATION_TPBORD_ABSTRACTGROUP_H

#include "LOCA_Continuation_AbstractGroup.H"    // base class

namespace LOCA { 
  
  namespace Bifurcation { 

    namespace TPBord {

      /*!  \brief Interface to underlying groups for turning point calculations
       */
      class AbstractGroup : 
	public virtual LOCA::Continuation::AbstractGroup {

      public:
  
	//! Default constructor.
	AbstractGroup() {}

	//! Destructor
	virtual ~AbstractGroup() {}
  
	//! Assignment operator
	virtual AbstractGroup& operator=(const AbstractGroup& source) = 0;

	//! applyJacobianInverse for nearly singular systems
	virtual NOX::Abstract::Group::ReturnType
	applySingularJacobianInverse(
				NOX::Parameter::List& params,
				const NOX::Abstract::Vector& input,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector& result) = 0;

	//! applyJacobianInverse for nearly singular systems, multiple RHS
	virtual NOX::Abstract::Group::ReturnType 
	applySingularJacobianInverseMulti(
			       NOX::Parameter::List& params,
			       const NOX::Abstract::Vector*const* inputs,
			       const NOX::Abstract::Vector& approxNullVec,
			       const NOX::Abstract::Vector& jacApproxNullVec,
			       NOX::Abstract::Vector** results,
			       int nVecs) = 0;

	//! Compute derivative of Jn w.r.t particular parameter param_id
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDp(const NOX::Abstract::Vector& nullVector,
		     const int param_id, NOX::Abstract::Vector& result) = 0;

	//! Compute derivative of Jn with respect to parameter, given JnVector
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDp(const NOX::Abstract::Vector& nullVector,
		     const int param_id,
		     const NOX::Abstract::Vector& JnVector,
		     NOX::Abstract::Vector& result) = 0;

	//! Compute vector (Jn)_{x}a  given vector a
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxa(const NOX::Abstract::Vector& nullVector,
		      const NOX::Abstract::Vector& aVector,
		      NOX::Abstract::Vector& result) = 0;
  
	//! Compute vector (Jn)_{x}a  given vector a, given JnVector
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxa(const NOX::Abstract::Vector& nullVector,
		      const NOX::Abstract::Vector& aVector,
		      const NOX::Abstract::Vector& JnVector,
		      NOX::Abstract::Vector& result) = 0;

	//! Compute vector (Jn)_{x}a+(Jn)_p given vector a, param_id 
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxaDp(const NOX::Abstract::Vector& nullVector,
			const NOX::Abstract::Vector& aVector,
			const int param_id,
			NOX::Abstract::Vector& result) = 0;

	//! Compute inner product between two vectors
	virtual double innerProduct(const NOX::Abstract::Vector& x,
				    const NOX::Abstract::Vector& y);
      };
    } // namespace TPBord
  } // namespace Continuation
} // namespace LOCA

#endif
