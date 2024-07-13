// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_TRANSPOSELINEARSYSTEM_ABSTRACTSTRATEGY_H
#define LOCA_EPETRA_TRANSPOSELINEARSYSTEM_ABSTRACTSTRATEGY_H

#include "Teuchos_RCP.hpp"

// Forward declares
class Epetra_Operator;
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  namespace Epetra {
    class Vector;
  }
}

namespace LOCA {

  namespace Epetra {

    //! Strategies for solving transpose linear systems
    namespace TransposeLinearSystem {

      //! A pure virtual interface for solving the transpose of a linear system
      class AbstractStrategy {

      public:

    //! Constructor
    AbstractStrategy() {}

    //! Destructor
    virtual ~AbstractStrategy() {}

    /*!
     * \brief Applies the inverse of the Jacobian matrix transpose to the
     * given input vector and puts the answer in result.
     *
     * Computes
     * \f[ v = J^{-T} u, \f]
     * where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector,
     * and \f$v\f$ is the result vector.
     *
     * The parameter list contains the linear solver options.
     */
    virtual bool
    applyJacobianTransposeInverse(Teuchos::ParameterList &params,
                      const NOX::Epetra::Vector &input,
                      NOX::Epetra::Vector &result) = 0;



    //! Evaluates the Jacobian-transpose based on the solution vector x.
    /*!
     * Note:  For flexibility, this method does not compute the original
     * Jacobian matrix.  It uses whatever is currently stored in the
     * linear system.
     */
    virtual bool
    createJacobianTranspose() = 0;

    /*!
     * \brief Explicitly constructs a preconditioner based on the solution
     * vector x and the parameter list p.
     */
    /*!
     * Note:  \c x is only needed for user-supplied preconditioners.
     * When using a built-in preconditioner (e.g., Ifpack), \c x will
     * note be used.
     */
    virtual bool
    createTransposePreconditioner(const NOX::Epetra::Vector& x,
                      Teuchos::ParameterList& p) = 0;

    //! Get Jacobian-transpose operator
    virtual Teuchos::RCP<Epetra_Operator>
    getJacobianTransposeOperator() = 0;

    //! Get transpose-preconditioner
    virtual Teuchos::RCP<Epetra_Operator>
    getTransposePreconditioner() = 0;

    //! Set Jacobian-transpose operator
    virtual void setJacobianTransposeOperator(
           const Teuchos::RCP<Epetra_Operator>& new_jac_trans) = 0;

    //! Set transpose-preconditioner
    virtual void setTransposePreconditioner(
          const Teuchos::RCP<Epetra_Operator>& new_prec_trans) = 0;

      }; // class AbstractStrategy

    } // namespace TransposeLinearSystem

  } // namespace Epetra

} // namespace LOCA


#endif
