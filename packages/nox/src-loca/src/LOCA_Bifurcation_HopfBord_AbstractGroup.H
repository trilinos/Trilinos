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

#ifndef LOCA_BIFURCATION_HOPFBORD_ABSTRACTGROUP_H
#define LOCA_BIFURCATION_HOPFBORD_ABSTRACTGROUP_H

#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"    // base class
#include "LOCA_TimeDependent_AbstractGroup.H"    // base class

namespace LOCA { 
  
  namespace Bifurcation { 

    namespace HopfBord {

      /*!  \brief Interface to underlying groups for Hopf calculations
       */
      class AbstractGroup : 
	public virtual LOCA::Bifurcation::TPBord::AbstractGroup,
	public virtual LOCA::TimeDependent::AbstractGroup {

      public:
  
	//! Default constructor.
	AbstractGroup() {}
 
	//! Destructor
	virtual ~AbstractGroup() {}

	//! Assignment operator
	virtual AbstractGroup& 
	operator=(const AbstractGroup& source) = 0;

	// Solve (J+iwB)*(y+iz) = (a+ib)
	virtual NOX::Abstract::Group::ReturnType
	applyComplexInverse(NOX::Parameter::List& params,
			    const NOX::Abstract::Vector& input_real,
			    const NOX::Abstract::Vector& input_imag,
			    double frequency,
			    NOX::Abstract::Vector& result_real,
			    NOX::Abstract::Vector& result_imag) const = 0;

	// Compute (J+iwB)*(y+iz)
	virtual NOX::Abstract::Group::ReturnType
	applyComplex(const NOX::Abstract::Vector& input_real,
		     const NOX::Abstract::Vector& input_imag,
		     double frequency,
		     NOX::Abstract::Vector& result_real,
		     NOX::Abstract::Vector& result_imag) const;

	//! Compute derivative of (J+iwB)(y+iz) with respect to particular 
	//! parameter param_id
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDp(const NOX::Abstract::Vector& yVector,
		     const NOX::Abstract::Vector& zVector,
		     double w,
		     const int param_id, 
		     NOX::Abstract::Vector& result_real,
		     NOX::Abstract::Vector& result_imag) = 0;

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
		     NOX::Abstract::Vector& result_imag) = 0;

	//! Compute vector (J+iwB)(y+iz))_{x}a  given vector a
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDxa(const NOX::Abstract::Vector& yVector,
		      const NOX::Abstract::Vector& zVector,
		      double w,
		      const NOX::Abstract::Vector& aVector,
		      NOX::Abstract::Vector& result_real,
		      NOX::Abstract::Vector& result_imag) = 0;
  
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
		      NOX::Abstract::Vector& result_imag) = 0;

	// Solve Complex system for multiple RHS
	virtual NOX::Abstract::Group::ReturnType
	applyComplexInverseMulti(
			     NOX::Parameter::List& params,
			     const NOX::Abstract::Vector* const* inputs_real,
			     const NOX::Abstract::Vector* const* inputs_imag,
			     double frequency,
			     NOX::Abstract::Vector** results_real,
			     NOX::Abstract::Vector** results_imag,
			     int nVecs) const;

      };
    } // namespace HopfBord
  } // namespace Continuation
} // namespace LOCA

#endif
