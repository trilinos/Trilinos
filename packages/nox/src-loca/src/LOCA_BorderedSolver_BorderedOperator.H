// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_BORDEREDSOLVER_BORDEREDOPERATOR_H
#define LOCA_BORDEREDSOLVER_BORDEREDOPERATOR_H

#include "Teuchos_RCP.hpp"

#include "LOCA_BorderedSolver_AbstractOperator.H" // base class
#include "LOCA_BorderedSystem_AbstractGroup.H"

namespace LOCA {

  namespace BorderedSolver {

    /*!
     * \brief Bordered solver operator representing as bordered Jacobian as
     * operator as implemented in the NOX::Abstract::Group.
     */
    class BorderedOperator : public LOCA::BorderedSolver::AbstractOperator {

    public:

      //! Constructor
      BorderedOperator(const Teuchos::RCP<const LOCA::BorderedSystem::AbstractGroup>& grp) : grpPtr(grp) {}

      //! Destructor
      virtual ~BorderedOperator() {}

      virtual Teuchos::RCP<const LOCA::BorderedSystem::AbstractGroup>
      getBorderedGroup() const { return grpPtr; }

      //! Apply the operator
      virtual NOX::Abstract::Group::ReturnType
      apply(const NOX::Abstract::MultiVector& /* X */,
        NOX::Abstract::MultiVector& /* Y */) const {
    return NOX::Abstract::Group::NotDefined;
      }

      //! Apply transpose of the operator
      virtual NOX::Abstract::Group::ReturnType
      applyTranspose(const NOX::Abstract::MultiVector& /* X */,
             NOX::Abstract::MultiVector& /* Y */) const {
    return NOX::Abstract::Group::NotDefined;
      }

      //! Apply inverse of the operator
      virtual NOX::Abstract::Group::ReturnType
      applyInverse(Teuchos::ParameterList& /* params */,
           const NOX::Abstract::MultiVector& /* B */,
           NOX::Abstract::MultiVector& /* X */) const {
    return NOX::Abstract::Group::NotDefined;
      }

      //! Apply inverse transpose of the operator
      virtual NOX::Abstract::Group::ReturnType
      applyInverseTranspose(Teuchos::ParameterList& /* params */,
                const NOX::Abstract::MultiVector& /* B */,
                NOX::Abstract::MultiVector& /* X */) const {
    return NOX::Abstract::Group::NotDefined;
      }

    private:

      //! Private to prohibit copying
      BorderedOperator(const BorderedOperator&);

      //! Private to prohibit copying
      BorderedOperator& operator = (const BorderedOperator&);

    protected:

      //! Pointer to group
      Teuchos::RCP<const LOCA::BorderedSystem::AbstractGroup> grpPtr;

    }; // Class BorderedOperator

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_BORDEREDOPERATOR_H
