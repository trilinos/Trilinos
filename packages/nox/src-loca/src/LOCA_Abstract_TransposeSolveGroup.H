// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_ABSTRACT_TRANSPOSESOLVEGROUP_H
#define LOCA_ABSTRACT_TRANSPOSESOLVEGROUP_H

#include "NOX_Abstract_Group.H"           // base class

namespace LOCA {

  namespace Abstract {

    /*!
     * \brief Abstract group interface class for solving the transpose
     * of the Jacobian.
     */
    /*!
     * This interface, derived from NOX::Abstract::Group, provides the
     * additional interface for solving the transpose of the Jacobian.
     */
    class TransposeSolveGroup : public virtual NOX::Abstract::Group {

    public:

      //! Constructor
      TransposeSolveGroup() {}

      //! Destructor
      virtual ~TransposeSolveGroup() {}

      //! Solve Jacobian-tranpose system
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianTransposeInverse(Teuchos::ParameterList& params,
                    const NOX::Abstract::Vector& input,
                    NOX::Abstract::Vector& result) const = 0;

      //! Solve Jacobian-tranpose system with multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianTransposeInverseMultiVector(
                Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector& input,
                NOX::Abstract::MultiVector& result) const = 0;

      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    }; // Class TransposeSolveGroup

  } // Namespace Abstract

} // Namespace LOCA

#endif // LOCA_ABSTRACT_TRANSPOSESOLVEGROUP_H
