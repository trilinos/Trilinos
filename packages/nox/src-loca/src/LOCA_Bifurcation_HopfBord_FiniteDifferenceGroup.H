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

#ifndef LOCA_BIFURCATION_HOPFBORD_FINITEDIFFERENCEGROUP_H
#define LOCA_BIFURCATION_HOPFBORD_FINITEDIFFERENCEGROUP_H

#include "LOCA_Bifurcation_HopfBord_AbstractGroup.H"          // base class
#include "LOCA_Bifurcation_TPBord_FiniteDifferenceGroup.H"    // base class

namespace LOCA { 
  
  namespace Bifurcation { 

    namespace HopfBord {

      /*!  \brief Interface to underlying groups for Hopf calculations
       */
      class FiniteDifferenceGroup : 
	public virtual LOCA::Bifurcation::HopfBord::AbstractGroup,
	public LOCA::Bifurcation::TPBord::FiniteDifferenceGroup {

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

	//! Compute derivative of (J+iwB)(y+iz) with respect to particular 
	//! parameter param_id
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDp(const NOX::Abstract::Vector& yVector,
		     const NOX::Abstract::Vector& zVector,
		     double w,
		     const int param_id, 
		     NOX::Abstract::Vector& result_real,
		     NOX::Abstract::Vector& result_imag);

	//! Compute derivative of (J+iwB)(y+iz) with respect to parameter, 
	//! given (J+iwB)(y+iz) vectors
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDp(const NOX::Abstract::Vector& yVector,
		     const NOX::Abstract::Vector& zVector,
		     double w,
		     const int param_id, 
		     const NOX::Abstract::Vector& Ce_real,
		     const NOX::Abstract::Vector& Ce_imag,
		     NOX::Abstract::Vector& result_real,
		     NOX::Abstract::Vector& result_imag);

	//! Compute vector (J+iwB)(y+iz))_{x}a  given vector a
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDxa(const NOX::Abstract::Vector& yVector,
		      const NOX::Abstract::Vector& zVector,
		      double w,
		      const NOX::Abstract::Vector& aVector,
		      NOX::Abstract::Vector& result_real,
		      NOX::Abstract::Vector& result_imag);
  
	//! Compute vector (J+iwB)(y+iz))_{x}a  given vector a, given 
	//! (J+iwB)(y+iz) vector
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDxa(const NOX::Abstract::Vector& yVector,
		      const NOX::Abstract::Vector& zVector,
		      double w,
		      const NOX::Abstract::Vector& aVector,
		      const NOX::Abstract::Vector& Ce_real,
		      const NOX::Abstract::Vector& Ce_imag,
		      NOX::Abstract::Vector& result_real,
		      NOX::Abstract::Vector& result_imag);
      };
    } // namespace HopfBord
  } // namespace Continuation
} // namespace LOCA

#endif
