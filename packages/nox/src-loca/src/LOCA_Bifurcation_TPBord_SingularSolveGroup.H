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

      /*!  
       * \brief Concrete class that provides concrete implementations of the
       * singular Jacobian solve methods of the 
       * LOCA::Bifurcation::TPBord::AbstractGroup.
       */
      /*!
       * The singular Jacobian solve methods are implemented by calling the
       * corresponding compute methods of the 
       * LOCA::SingularJacobianSolve::Manager class.  It is through this class,
       * and the passed parameter list in the constructor, that the particular
       * algorithm for performing the singular solves is chosen.  See the
       * LOCA::SingularJacobianSolve::Manager class description for a 
       * description of the relevant parameters and the corresponding
       * singular solve algorithms.  In the future, the functionality provided
       * by the manager class will be directly incorporated into this class.
       */
      class SingularSolveGroup : 
	public virtual LOCA::Bifurcation::TPBord::AbstractGroup {

      public:

	//! Constructor with const parameter list
	SingularSolveGroup(const NOX::Parameter::List& params = NOX::Parameter::List());
  
	//! Constructor with non-const parameter list
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
	/*!
	 * This method is to compute the solution to \f$Jx=b\f$ for nearly
	 * singular \f$J\f$.  \b approxNullVec contains an approximation to the
	 * null vector \f$n\f$ such that \f$Jn\approx0\f$, and 
	 * \b jacApproxNullVec
	 * contains \f$Jn\f$.  The idea is these quantities are supposed to
	 * help in obtaining an accurate solution to \f$Jx=b\f$.
	 *
	 * The method is implemented by calling the 
	 * LOCA::SingularJacobianSolve::Manager::compute() method of the 
	 * LOCA::SingularJacobianSolve::Manager class.
	 */
	virtual NOX::Abstract::Group::ReturnType
	applySingularJacobianInverse(
				NOX::Parameter::List& params,
				const NOX::Abstract::Vector& input,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector& result);

	//! Version of applySingularJacobianInverse() for multiple RHS.
	/*!
	 * The method is implemented by calling the 
	 * LOCA::SingularJacobianSolve::Manager::computeMulti() method of the 
	 * LOCA::SingularJacobianSolve::Manager class.
	 */
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
  } // namespace Bifurcation
} // namespace LOCA

#endif
