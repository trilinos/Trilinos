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

#ifndef LOCA_BIFURCATION_TPBORD_SINGULARSOLVEGROUP_H
#define LOCA_BIFURCATION_TPBORD_SINGULARSOLVEGROUP_H

#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"    // base class
#include "LOCA_SingularJacobianSolve_Manager.H"       // class data element
#include "NOX_Parameter_List.H"

namespace LOCA { 
  
  namespace Bifurcation { 

    namespace TPBord {

      /*!  \brief Interface to underlying groups for turning point calculations
       */
      class SingularSolveGroup : 
	public virtual LOCA::Bifurcation::TPBord::AbstractGroup {

      public:

	//! Default constructor.
	SingularSolveGroup(const NOX::Parameter::List& params = NOX::Parameter::List());
  
	//! Default constructor.
	SingularSolveGroup(NOX::Parameter::List& params);
 
	//! Copy constructor
	SingularSolveGroup(const SingularSolveGroup& source, 
			   NOX::CopyType type = NOX::DeepCopy);

	//! Destructor
	virtual ~SingularSolveGroup();
  
	//! Assignment operator
	virtual SingularSolveGroup& 
	operator=(const SingularSolveGroup& source);

	//! applyJacobianInverse for nearly singular systems
	virtual NOX::Abstract::Group::ReturnType
	applySingularJacobianInverse(
				NOX::Parameter::List& params,
				const NOX::Abstract::Vector& input,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector& result);

	//! applyJacobianInverse for nearly singular systems, multiple RHS
	virtual NOX::Abstract::Group::ReturnType 
	applySingularJacobianInverseMulti(
			       NOX::Parameter::List& params,
			       const NOX::Abstract::Vector*const* inputs,
			       const NOX::Abstract::Vector& approxNullVec,
			       const NOX::Abstract::Vector& jacApproxNullVec,
			       NOX::Abstract::Vector** results,
			       int nVecs);

      protected:

	//! Manager for applySingularJacobianInverse[Multi]
	LOCA::SingularJacobianSolve::Manager singularSolveManager;
      };
    } // namespace TPBord
  } // namespace Continuation
} // namespace LOCA

#endif
