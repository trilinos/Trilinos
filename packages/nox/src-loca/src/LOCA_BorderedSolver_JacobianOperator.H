// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_BORDEREDSOLVER_JACOBIANOPERATOR_H
#define LOCA_BORDEREDSOLVER_JACOBIANOPERATOR_H

#include "Teuchos_RCP.hpp"

#include "LOCA_BorderedSolver_AbstractOperator.H" // base class

namespace LOCA {

  namespace BorderedSolver {

    /*!
     * \brief Bordered solver operator representing the Jacobian as
     * implemented in the NOX::Abstract::Group.
     */
    class JacobianOperator : public LOCA::BorderedSolver::AbstractOperator {

    public:

      //! Constructor
      JacobianOperator(
         const Teuchos::RCP<const NOX::Abstract::Group>& grp);

      //! Destructor
      virtual ~JacobianOperator();

      //! Get group pointer
      virtual Teuchos::RCP<const NOX::Abstract::Group>
      getGroup() const;

      //! Apply the operator
      virtual NOX::Abstract::Group::ReturnType
      apply(const NOX::Abstract::MultiVector& X,
        NOX::Abstract::MultiVector& Y) const;

      //! Apply transpose of the operator
      virtual NOX::Abstract::Group::ReturnType
      applyTranspose(const NOX::Abstract::MultiVector& X,
             NOX::Abstract::MultiVector& Y) const;

      //! Apply inverse of the operator
      virtual NOX::Abstract::Group::ReturnType
      applyInverse(Teuchos::ParameterList& params,
           const NOX::Abstract::MultiVector& B,
           NOX::Abstract::MultiVector& X) const;

      //! Apply inverse transpose of the operator
      /*!
       * Group must be of type LOCA::Abstract::TransposeSolveGroup
       * for this method to be defined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyInverseTranspose(Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector& B,
                NOX::Abstract::MultiVector& X) const;

    private:

      //! Private to prohibit copying
      JacobianOperator(const JacobianOperator&);

      //! Private to prohibit copying
      JacobianOperator& operator = (const JacobianOperator&);

    protected:

      //! Pointer to group
      Teuchos::RCP<const NOX::Abstract::Group> grpPtr;

    }; // Class JacobianOperator

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_JACOBIANOPERATOR_H
