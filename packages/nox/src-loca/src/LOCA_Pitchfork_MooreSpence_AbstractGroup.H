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

#ifndef LOCA_PITCHFORK_MOORESPENCE_ABSTRACTGROUP_H
#define LOCA_PITCHFORK_MOORESPENCE_ABSTRACTGROUP_H

#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"    // base class

namespace LOCA { 
  
  //! Groups and vectors for pitchfork bifurcations
  namespace Pitchfork { 

    /*! 
     * \brief Groups and vectors for locating pitchfork bifurcations
     * using the Moore-Spence formulation.
     */
    namespace MooreSpence {

      /*! 
       * \brief Interface to underlying groups for pitchfork calculations
       * using the Moore-Spence formulation.
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate pitchforks using the bordering algorithms for the 
       * Moore-Spence pitchfork formulation (see 
       * LOCA::Pitchfork::MooreSpence::ExtendedGroup for a description of the 
       * governing equations).  
       *
       * This class is derived from the 
       * LOCA::TurningPoint::MooreSpence::AbstractGroup and declares a single
       * virtual method, innerProduct(), to compute the inner product of a 
       * vector with the asymmetry vector.  It has a default implementation
       * given by the dot product, but should be overloaded for any 
       * problem that has a different definition for the inner product.
       */
      class AbstractGroup : 
	public virtual LOCA::TurningPoint::MooreSpence::AbstractGroup {

      public:
  
	//! Default constructor.
	AbstractGroup() {}

	//! Destructor
	virtual ~AbstractGroup() {}

	//! Compute the inner product of \c a and \c b.
	/*!
	 * The default implementation is given by the dot product of
	 * \c a and \c b.
	 */
	virtual double innerProduct(const NOX::Abstract::Vector& a,
				    const NOX::Abstract::Vector& b) const {
	  return a.innerProduct(b);
	}

	//! Compute the inner product of \c a and \c b.
	/*!
	 * The default implementation is given by the dot product of
	 * \c a and \c b.
	 */
	virtual void innerProduct(
			    const NOX::Abstract::MultiVector& a,
			    const NOX::Abstract::MultiVector& b,
			    NOX::Abstract::MultiVector::DenseMatrix& c) const {
	  b.multiply(1.0, a, c);
	}

	//! Bring NOX::Abstract::Group::operator=() into scope
	using NOX::Abstract::Group::operator=;

      }; // class AbstractGroup

    } // namespace MooreSpence

  } // namespace Pitchfork

} // namespace LOCA

#endif
