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
     * Hopf calculations).  It provides pure virtual methods for computing
     * and manipulating the shifted matrix \f$\alpha J + \beta M\f$ where
     * \f$J\f$ is the Jacobian matrix and \f$M\f$ is the mass matrix.
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

      //! Compute the shifted matrix
      virtual NOX::Abstract::Group::ReturnType
      computeShiftedMatrix(double alpha, double beta) = 0;

      //! Multiply the shifted matrix by a vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrix(const NOX::Abstract::Vector& input,
			 NOX::Abstract::Vector& result) const = 0;

      //! Multiply the shifted matrix by a multi-vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixMultiVector(
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const = 0;
      
      /*! 
       * \brief Apply the inverse of the shifted matrix by a multi-vector, as 
       * needed by the shift-and-invert and generalized Cayley transformations.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixInverseMultiVector(
			        Teuchos::ParameterList& params, 
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const = 0;


      //@}

      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    };
  } // namespace TimeDependent
} // namespace LOCA

#endif
