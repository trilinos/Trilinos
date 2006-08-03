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

#ifndef LOCA_TURNINGPOINT_MOORESPENCE_FINITEDIFFERENCEGROUP_H
#define LOCA_TURNINGPOINT_MOORESPENCE_FINITEDIFFERENCEGROUP_H

#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"   // base class
#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H" // base class

namespace LOCA { 
  
  namespace TurningPoint { 

    namespace MooreSpence {

      /*!  
       * \brief Concrete class that provides concrete implementations of the
       * derivative computation methods of the 
       * LOCA::TurningPoint::MooreSpence::AbstractGroup using first-order 
       * finite differencing.
       */
      /*! 
       * The finite-differencing calculations are actually implemented by
       * the LOCA::DerivUtils class, and a custom DerivUtils object can be
       * passed through the constructor of this class.  However, in the future
       * the calculations encapsulated in the DerivUtils class may be 
       * incorporated directly into this class and other finite-differencing
       * child classes.
       */
      class FiniteDifferenceGroup : 
	public virtual LOCA::TurningPoint::MooreSpence::AbstractGroup,
	public virtual LOCA::MultiContinuation::FiniteDifferenceGroup {

      public:
  
	//! Constructor.
	FiniteDifferenceGroup();
 
	//! Copy constructor
	FiniteDifferenceGroup(const FiniteDifferenceGroup& source, 
			      NOX::CopyType type = NOX::DeepCopy);

	//! Destructor
	virtual ~FiniteDifferenceGroup();

	/*! 
	 * \brief Computes the derivative \f$\partial Jn/\partial p\f$ 
	 */
	/*! 
	 * The calculation is implemented by calling the corresponding
	 * LOCA::DerivUtils::computeDJnDp() method of the passed 
	 * LOCA::DerivUtils object.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDpMulti(const vector<int>& paramIDs, 
			  const NOX::Abstract::Vector& nullVector,
			  NOX::Abstract::MultiVector& result,
			  bool isValid);

	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction 
	 * \f$a\f$. 
	 */
	/*! 
	 * The calculation is implemented by calling the corresponding
	 * LOCA::DerivUtils::computeDJnDxa() method of the passed 
	 * LOCA::DerivUtils object.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxaMulti(const NOX::Abstract::Vector& nullVector,
			   const NOX::Abstract::MultiVector& aVector,
			   NOX::Abstract::MultiVector& result);

	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction 
	 * \f$a\f$. 
	 */
	/*! 
	 * The calculation is implemented by calling the corresponding
	 * LOCA::DerivUtils::computeDJnDxa() method of the passed 
	 * LOCA::DerivUtils object.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxaMulti(const NOX::Abstract::Vector& nullVector,
			   const NOX::Abstract::Vector& JnVector,
			   const NOX::Abstract::MultiVector& aVector,
			   NOX::Abstract::MultiVector& result);

	//! Bring NOX::Abstract::Group::operator=() into scope
	using NOX::Abstract::Group::operator=;
  
      }; // class FiniteDifferenceGroup
    } // namespace MooreSpence
  } // namespace TurningPoint
} // namespace LOCA

#endif
