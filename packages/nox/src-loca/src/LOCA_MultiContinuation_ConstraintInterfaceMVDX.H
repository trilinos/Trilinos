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

#ifndef LOCA_MULTICONTINUATION_CONSTRAINTINTERFACEMVDX_H
#define LOCA_MULTICONTINUATION_CONSTRAINTINTERFACEMVDX_H

#include "LOCA_MultiContinuation_ConstraintInterface.H" // base class

namespace LOCA {

  namespace MultiContinuation {

    /*! 
     * \brief Abstract interface for the constraint portion of a constrained 
     * nonlinear system for constraints that support computing a solution
     * component derivative as a multi-vector.
     */
    /*!
     * This class extends the LOCA::MultiContinuation::ConstraintInterface
     * to support constraints that support computing the entire derivative
     * with respect to the solution components (x) and storing the resulting
     * derivative as a multivector.  This interface adds one additional
     * method, getConstraintDerivativesX(), that returns this derivative.
     * Additionally, it implements the applyConstraintDerivativesX() methods
     * using standard multi-vector operations.
     */
    class ConstraintInterfaceMVDX : 
      public virtual LOCA::MultiContinuation::ConstraintInterface {

    public:

      //! Constructor
      ConstraintInterfaceMVDX() {}

      //! Destructor
      virtual ~ConstraintInterfaceMVDX() {}

      //! Return solution component of constraint derivatives
      /*!
       * May return NULL if constraint derivative is zero
       */
      virtual const NOX::Abstract::MultiVector*
      getDX() const = 0;

      //! Compute result_p = alpha * dg/dx * input_x
      /*!
       * This method is implemented using getConstraintDerivativesX()
       * and the NOX::Abstract::MultiVector::multiply() method.
       */
      virtual NOX::Abstract::Group::ReturnType
      multiplyDX(double alpha, 
		 const NOX::Abstract::MultiVector& input_x,
		 NOX::Abstract::MultiVector::DenseMatrix& result_p) const;

      //! Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x
      /*!
       * This method is implemented using getConstraintDerivativesX()
       * and the NOX::Abstract::MultiVector::update() method.
       */
      virtual NOX::Abstract::Group::ReturnType
      addDX(Teuchos::ETransp transb,
	    double alpha, 
	    const NOX::Abstract::MultiVector::DenseMatrix& b,
	    double beta,
	    NOX::Abstract::MultiVector& result_x) const;

    }; // Class ConstraintInterfaceMVDX

  } // Namespace MultiContinuation

} // Namespace LOCA

#endif // LOCA_MULTICONTINUATION_CONSTRAINTINTERFACE_H
