// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_BORDEREDSOLVER_ABSTRACTOPERATOR_H
#define LOCA_BORDEREDSOLVER_ABSTRACTOPERATOR_H

#include "NOX_Abstract_Group.H"           // for ReturnType
#include "NOX_Abstract_MultiVector.H"

namespace LOCA {

  namespace BorderedSolver {

    /*!
     * \brief Abstract interface class representing an operator for solving
     * bordered sets of linear equations */
    class AbstractOperator {

    public:

      //! Constructor
      AbstractOperator() {}

      //! Destructor
      virtual ~AbstractOperator() {}

      //! Apply the operator
      virtual NOX::Abstract::Group::ReturnType
      apply(const NOX::Abstract::MultiVector& X,
        NOX::Abstract::MultiVector& Y) const = 0;

      //! Apply transpose of the operator
      virtual NOX::Abstract::Group::ReturnType
      applyTranspose(const NOX::Abstract::MultiVector& X,
             NOX::Abstract::MultiVector& Y) const = 0;

      //! Apply inverse of the operator
      virtual NOX::Abstract::Group::ReturnType
      applyInverse(Teuchos::ParameterList& params,
           const NOX::Abstract::MultiVector& B,
           NOX::Abstract::MultiVector& X) const = 0;

      //! Apply inverse transpose of the operator
      virtual NOX::Abstract::Group::ReturnType
      applyInverseTranspose(Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector& B,
                NOX::Abstract::MultiVector& X) const = 0;

    private:

      //! Private to prohibit copying
      AbstractOperator(const AbstractOperator&);

      //! Private to prohibit copying
      AbstractOperator& operator = (const AbstractOperator&);

    }; // Class AbstractOperator

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_ABSTRACTOPERATOR_H
