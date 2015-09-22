// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
