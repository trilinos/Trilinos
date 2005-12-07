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

#ifndef LOCA_TURNINGPOINT_MOORESPENCE_ABSTRACTGROUP_H
#define LOCA_TURNINGPOINT_MOORESPENCE_ABSTRACTGROUP_H

#include "LOCA_MultiContinuation_AbstractGroup.H"    // base class

namespace LOCA { 
  
  //! Groups and vectors for turning point bifurcations
  namespace TurningPoint { 

    /*! 
     * \brief Groups and vectors for locating turning point bifurcations
     * using the Moore-Spence formulation.
     */
    namespace MooreSpence {

      /*! 
       * \brief Interface to underlying groups for turning point calculations
       * using the Moore-Spence formulation.
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate turning points and pitchforks using the bordering
       * algorithms for the Moore-Spence turning point and pitchfork 
       * formulations (see LOCA::TurningPoint::MooreSpence::ExtendedGroup
       * for a description of the governing equations).  
       *
       * This class
       * is derived from the LOCA::Continuation::AbstractGroup and declares
       * several pure virtual methods compute various 
       * derivatives of \f$Jn\f$ for a given \f$n\f$.  Default implementations
       * for the derivativesusing finite differencing are implemented in the 
       * LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup.
       *
       * The last virtual method, innerProduct() has a default implementation
       * given by the dot product, but should be overloaded for any 
       * problem that has a different definition for the inner product (this
       * is only used in the LOCA::Bifurcation::PitchforkBord::ExtendedGroup
       * for the asymmetry).
       */
      class AbstractGroup : 
	public virtual LOCA::MultiContinuation::AbstractGroup {

      public:
  
	//! Default constructor.
	AbstractGroup() {}

	//! Destructor
	virtual ~AbstractGroup() {}

	/*! 
	 * @name Pure virtual methods 
	 * These methods must be defined by any concrete implementation
	 */
	//@{

	/*! 
	 * \brief Computes the derivative \f$\partial Jn/\partial p\f$ 
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDpMulti(const vector<int>& paramIDs, 
			  const NOX::Abstract::Vector& nullVector,
			  NOX::Abstract::MultiVector& result,
			  bool isValid) = 0;

	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction 
	 * \f$a\f$. 
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxaMulti(const NOX::Abstract::Vector& nullVector,
			   const NOX::Abstract::MultiVector& aVector,
			   NOX::Abstract::MultiVector& result) = 0;

	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction 
	 * \f$a\f$. 
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxaMulti(const NOX::Abstract::Vector& nullVector,
			   const NOX::Abstract::Vector& JnVector,
			   const NOX::Abstract::MultiVector& aVector,
			   NOX::Abstract::MultiVector& result) = 0;

	//@}

	//! Bring NOX::Abstract::Group::operator=() into scope
	using NOX::Abstract::Group::operator=;

      }; // class AbstractGroup

    } // namespace MooreSpence

  } // namespace TurningPoint

} // namespace LOCA

#endif
