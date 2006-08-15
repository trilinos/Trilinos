// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_HOPF_MOORESPENCE_FINITEDIFFERENCEGROUP_H
#define LOCA_HOPF_MOORESPENCE_FINITEDIFFERENCEGROUP_H

#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"                 // base class
#include "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H" // base class

namespace LOCA { 
  
  namespace Hopf { 

    namespace MooreSpence {

      /*!  
       * \brief Concrete class that provides concrete implementations of the
       * derivative computation methods of the 
       * LOCA::Hopf::MooreSpence::AbstractGroup using first-order finite 
       * differencing.
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
	public virtual LOCA::Hopf::MooreSpence::AbstractGroup,
	public virtual LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup {

      public:
  
	//! Constructor.
	FiniteDifferenceGroup();
 
	//! Copy constructor
	FiniteDifferenceGroup(const FiniteDifferenceGroup& source, 
			      NOX::CopyType type = NOX::DeepCopy);

	//! Destructor
	virtual ~FiniteDifferenceGroup();

	/*! 
	 * \brief Computes the derivative 
	 * \f$\frac{\partial (J+i\omega B)(y+iz)}{\partial p}\f$ 
	 * where \f$p\f$ is the parameter indexed by \b paramIDs.
	 */
	/*! 
	 * The calculation is implemented by calling the corresponding
	 * LOCA::DerivUtils::computeDCeDp() method of the passed 
	 * LOCA::DerivUtils object.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDp(const vector<int>& paramIDs, 
		     const NOX::Abstract::Vector& yVector,
		     const NOX::Abstract::Vector& zVector,
		     double w,
		     NOX::Abstract::MultiVector& result_real,
		     NOX::Abstract::MultiVector& result_imag,
		     bool isValid);

	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial (J+i\omega B)(y+iz)}{\partial x} a\f$ for the 
	 * given direction \f$a\f$. 
	 */
	/*! 
	 * The calculation is implemented by calling the corresponding
	 * LOCA::DerivUtils::computeDCeDxa() method of the passed 
	 * LOCA::DerivUtils object.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDxa(const NOX::Abstract::Vector& yVector,
		      const NOX::Abstract::Vector& zVector,
		      double w,
		      const NOX::Abstract::MultiVector& aVector,
		      NOX::Abstract::MultiVector& result_real,
		      NOX::Abstract::MultiVector& result_imag);
  
	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial (J+i\omega B)(y+iz)}{\partial x} a\f$ for the 
	 * given direction \f$a\f$. The 
	 * arguments \b Ce_real and \b Ce_imag hold the real and imaginary
	 * components of \f$(J+i\omega B)(y+iz)\f$.
	 */
	/*! 
	 * The calculation is implemented by calling the corresponding
	 * LOCA::DerivUtils::computeDCeDxa() method of the passed 
	 * LOCA::DerivUtils object.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDxa(const NOX::Abstract::Vector& yVector,
		      const NOX::Abstract::Vector& zVector,
		      double w,
		      const NOX::Abstract::MultiVector& aVector,
		      const NOX::Abstract::Vector& Ce_real,
		      const NOX::Abstract::Vector& Ce_imag,
		      NOX::Abstract::MultiVector& result_real,
		      NOX::Abstract::MultiVector& result_imag);

	//! Bring NOX::Abstract::Group::operator=() into scope
	using NOX::Abstract::Group::operator=;

      };
    } // namespace HopfBord
  } // namespace Bifurcation
} // namespace LOCA

#endif
