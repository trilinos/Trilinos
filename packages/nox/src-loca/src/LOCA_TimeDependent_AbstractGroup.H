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

#ifndef LOCA_TIMEDEPENDENT_ABSTRACTGROUP_H
#define LOCA_TIMEDEPENDENT_ABSTRACTGROUP_H

#include "LOCA_MultiContinuation_AbstractGroup.H"    // base class

namespace LOCA { 

  namespace TimeDependent {

    /*!  
     * \brief Interface to underlying groups for time dependent systems
     */
    /*!
     * This abstract class provides an interface for time dependent problems,
     * i.e., problems with a mass matrix (typically used in eignvalue or
     * Hopf calculations).  It provides pure virtual methods for manipulating
     * the mass matrix and adding the mass and Jacobian matrices.
     *
     * Where possible, default implementations are provided assuming the
     * mass matrix is the identity.
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
       * \brief Apply the inverse of the shifted matrix by a vector, as 
       * needed by the shift-and-invert and generalized Cayley transformations.
       */
      /*!
       * Note that because of the variety of ways this method can be 
       * implemented, the interface for this method is different than 
       * other applyInverse() methods.  Namely, the method is not const,
       * and there are no computeShiftMatrix(), isShiftedMatrix() methods.
       * Rather, it is assumed the group will compute the shifted matrix
       * if necessary, and then apply the solve.  It is up to the 
       * implementing group to determine if the shifted matrix needs to be
       * recomputed.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixInverse(Teuchos::ParameterList& params, 
				const NOX::Abstract::Vector& input,
                                NOX::Abstract::Vector& result,
                                double shift) = 0;


      //@}

      /*! 
       * @name Virtual methods with default implementations 
       * These methods should be overloaded in a concrete implementation
       * if more appropriate/efficient approaches are available.
       */
      //@{

      //! Compute the mass matrix of the system
      /*!
       * The default assumes the mass matrix is the identity, and 
       * therefore merely returns NOX::Abstract::Group::Ok.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeMassMatrix();

      //! Is mass matrix valid
      /*!
       * The default assumes the mass matrix is the identity, and 
       * therefore merely returns true.
       */
      virtual bool isMassMatrix() const;

      //! Multiply the mass matrix by a vector.
      /*!
       * The default assumes the mass matrix is the identity, and 
       * therefore merely sets \em result = \em input.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyMassMatrix(const NOX::Abstract::Vector& input,
		      NOX::Abstract::Vector& result) const;

      //! Multiply the shifted matrix by a vector
      /*!
       * The default assumes the mass matrix is the identity, and 
       * therefore merely sets 
       * \em result = \em J * \em input + \em shift * \em input.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrix(const NOX::Abstract::Vector& input,
                         NOX::Abstract::Vector& result,
                         double shift) const;

      //@}

      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    };
  } // namespace TimeDependent
} // namespace LOCA

#endif
