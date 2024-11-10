// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_MULTICONTINUATION_COMPOSITECONSTRAINTMVDX_H
#define LOCA_MULTICONTINUATION_COMPOSITECONSTRAINTMVDX_H

#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H" // base class
#include "LOCA_MultiContinuation_CompositeConstraint.H"     // base class

namespace LOCA {

  namespace MultiContinuation {

    /*!
     * \brief Implementation of
     * LOCA::MultiContinuation::ConstraintInterfaceMVDX for composite
     * constraints, i.e., a constraint comprised of multiple, separate
     * constraints.
     */
    class CompositeConstraintMVDX :
      public virtual LOCA::MultiContinuation::ConstraintInterfaceMVDX,
      public virtual LOCA::MultiContinuation::CompositeConstraint {

    public:

      //! Constructor
      CompositeConstraintMVDX(
       const Teuchos::RCP<LOCA::GlobalData>& global_data,
       const std::vector< Teuchos::RCP<
       LOCA::MultiContinuation::ConstraintInterfaceMVDX> >& constraintObjects);

      //! Copy constructor
      CompositeConstraintMVDX(const CompositeConstraintMVDX& source,
                  NOX::CopyType type = NOX::DeepCopy);

      //! Destructor
      ~CompositeConstraintMVDX();

      /*!
       * @name Implementation of LOCA::MultiContinuation::ConstraintInterface
       * virtual methods
       */
      //@{

      //! Copy
      virtual void copy(const ConstraintInterface& source);

      //! Cloning function
      virtual
      Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Compute derivative of constraints w.r.t. solution vector x
      virtual NOX::Abstract::Group::ReturnType
      computeDX();

      //! Compute result_p = alpha * dg/dx * input_x
      /*!
       * Note that if there are n constraints and input_x has m columns,
       * result_p should be a n by m matrix and is equivalent to
       * \code
       *    input_x.multiply(alpha, dgdx, result_p).
       * \endcode
       */
      virtual NOX::Abstract::Group::ReturnType
      multiplyDX(double alpha,
          const NOX::Abstract::MultiVector& input_x,
          NOX::Abstract::MultiVector::DenseMatrix& result_p) const;

      //! Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x
      /*!
       * Note that this should be equivalent to
       * \code
       *    result_x.update(transb, alpha, dgdx, b, beta);
       * \endcode
       */
      virtual NOX::Abstract::Group::ReturnType
      addDX(Teuchos::ETransp transb,
        double alpha,
        const NOX::Abstract::MultiVector::DenseMatrix& b,
        double beta,
        NOX::Abstract::MultiVector& result_x) const;

      //@}

      /*!
       * @name Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX
       * virtual methods
       */
      //@{

      //! Return solution component of constraint derivatives
      virtual const NOX::Abstract::MultiVector*
      getDX() const;

      //@}

    private:

      //! Prohibit generation and use of operator=()
      CompositeConstraintMVDX&
      operator=(const CompositeConstraintMVDX& source);

    protected:

      //! Array of constraintMVDX pointers
      std::vector< Teuchos::RCP<
    LOCA::MultiContinuation::ConstraintInterfaceMVDX> > constraintMVDXPtrs;

      //! Composite constraint derivative
      Teuchos::RCP<NOX::Abstract::MultiVector> compositeDX;

    }; // Class CompositeConstraintMVDX

  } // namespace MultiContinuation

} // namespace LOCA

#endif // LOCA_MULTICONTINUATION_COMPOSITECONSTRAINT_H
