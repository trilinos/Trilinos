// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
