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

#ifndef LOCA_BIFURCATION_TPBORD_FINITEDIFFERENCEGROUP_H
#define LOCA_BIFURCATION_TPBORD_FINITEDIFFERENCEGROUP_H

#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"   // base class
#include "LOCA_Continuation_FiniteDifferenceGroup.H" // base class

namespace LOCA { 
  
  namespace Bifurcation { 

    namespace TPBord {

      /*!  \brief Interface to underlying groups for turning point calculations
       */
      class FiniteDifferenceGroup : 
	public virtual LOCA::Bifurcation::TPBord::AbstractGroup,
	public LOCA::Continuation::FiniteDifferenceGroup {

      public:
  
	//! Default constructor.
	FiniteDifferenceGroup(const LOCA::DerivUtils& deriv);
 
	//! Copy constructor
	FiniteDifferenceGroup(const FiniteDifferenceGroup& source, 
			      NOX::CopyType type = NOX::DeepCopy);

	//! Destructor
	virtual ~FiniteDifferenceGroup();
  
	//! Assignment operator
	virtual FiniteDifferenceGroup& 
	operator=(const FiniteDifferenceGroup& source);

	//! Compute derivative of Jn w.r.t particular parameter param_id
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDp(const NOX::Abstract::Vector& nullVector,
		     const int param_id, NOX::Abstract::Vector& result);

	//! Compute derivative of Jn with respect to parameter, given JnVector
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDp(const NOX::Abstract::Vector& nullVector,
		     const int param_id,
		     const NOX::Abstract::Vector& JnVector,
		     NOX::Abstract::Vector& result);

	//! Compute vector (Jn)_{x}a  given vector a
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxa(const NOX::Abstract::Vector& nullVector,
		      const NOX::Abstract::Vector& aVector,
		      NOX::Abstract::Vector& result);
  
	//! Compute vector (Jn)_{x}a  given vector a, given JnVector
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxa(const NOX::Abstract::Vector& nullVector,
		      const NOX::Abstract::Vector& aVector,
		      const NOX::Abstract::Vector& JnVector,
		      NOX::Abstract::Vector& result);

	//! Compute vector (Jn)_{x}a+(Jn)_p given vector a, param_id 
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxaDp(const NOX::Abstract::Vector& nullVector,
			const NOX::Abstract::Vector& aVector,
			const int param_id,
			NOX::Abstract::Vector& result);
      };
    } // namespace TPBord
  } // namespace Continuation
} // namespace LOCA

#endif
