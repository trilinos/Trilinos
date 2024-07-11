// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_TRANSPOSELINEARSYSTEM_TRANSPOSEPRECONDITIONER_H
#define LOCA_EPETRA_TRANSPOSELINEARSYSTEM_TRANSPOSEPRECONDITIONER_H

#include "LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H" // base class

// forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  namespace Epetra {
    class LinearSystem;
    class Scaling;
  }
}
namespace LOCA {
  class GlobalData;
}

namespace LOCA {

  namespace Epetra {

    namespace TransposeLinearSystem {

      /*!
       * \brief Method for solving the transpose of a linear system by
       * using the transpose of the preconditioner.
       */
      /* The class is intialized via the \c solverParams parameter list
       * argument to the constructor.  The parameters this class recognizes
       * are:
       * <ul>
       * <li> "Transpose Scaling"  --
       *    [Teuchos::RCP<NOX::Epetra::Scaling>]  -
       *    Scaling object for transpose solve
       * </ul>
       */
      class TransposePreconditioner : public virtual AbstractStrategy {

      public:

    //! Constructor
    TransposePreconditioner(
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
        const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
        const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys);

    //! Destructor
    virtual ~TransposePreconditioner();

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
                      NOX::Epetra::Vector &result);



    //! Evaluates the Jacobian-transpose based on the solution vector x.
    /*!
     * Note:  For flexibility, this method does not compute the original
     * Jacobian matrix.  It uses whatever is currently stored in the
     * linear system.
     */
    virtual bool
    createJacobianTranspose();

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
                      Teuchos::ParameterList& p);

    //! Get Jacobian-transpose operator
    virtual Teuchos::RCP<Epetra_Operator>
    getJacobianTransposeOperator();

    //! Get transpose-preconditioner
    virtual Teuchos::RCP<Epetra_Operator>
    getTransposePreconditioner();

    //! Set Jacobian-transpose operator
    virtual void setJacobianTransposeOperator(
           const Teuchos::RCP<Epetra_Operator>& new_jac_trans);

    //! Set transpose-preconditioner
    virtual void setTransposePreconditioner(
          const Teuchos::RCP<Epetra_Operator>& new_prec_trans);

      protected:

    //! Global data object
    Teuchos::RCP<LOCA::GlobalData> globalData;

    //! Linear system solver
    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;

    //! Jacobian operator
    Teuchos::RCP<Epetra_Operator> jac;

    //! Preconditioner
    Teuchos::RCP<Epetra_Operator> prec;

    //! Scaling object for transpose solve
    Teuchos::RCP<NOX::Epetra::Scaling> scaling_trans;

      }; // class TransposePreconditioner

    } // namespace TransposeLinearSystem

  } // namespace Epetra

} // namespace LOCA


#endif
