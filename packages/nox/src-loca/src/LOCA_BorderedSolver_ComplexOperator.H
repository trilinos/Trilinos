// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_BORDEREDSOLVER_COMPLEXOPERATOR_H
#define LOCA_BORDEREDSOLVER_COMPLEXOPERATOR_H

#include "Teuchos_RCP.hpp"

#include "LOCA_BorderedSolver_AbstractOperator.H" // base class

// forward declarations
namespace LOCA {
  namespace Hopf {
    namespace MooreSpence {
      class AbstractGroup;
    }
  }
}

namespace LOCA {

  namespace BorderedSolver {

    /*!
     * \brief Bordered solver operator representing the \f$J + i\omega M\f$ as
     * implemented in the LOCA::Hopf::MooreSpence::AbstractGroup.
     */
    class ComplexOperator : public LOCA::BorderedSolver::AbstractOperator {

    public:

      //! Constructor
      ComplexOperator(const Teuchos::RCP<const LOCA::Hopf::MooreSpence::AbstractGroup>& grp,
              double omega);

      //! Destructor
      virtual ~ComplexOperator();

      //! Get group pointer
      virtual Teuchos::RCP<const NOX::Abstract::Group>
      getGroup() const;

      //! Get frequency
      virtual double getFrequency() const;

      //! Apply the operator
      virtual NOX::Abstract::Group::ReturnType
      apply(const NOX::Abstract::MultiVector& X,
        NOX::Abstract::MultiVector& Y) const;

      //! Apply transpose of the operator
      /*!
       * Group must be of type LOCA::Hopf::MinimallyAugmented::AbstractGroup
       * for this method to be implemented.
       */
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
       * Group must be of type LOCA::Hopf::MinimallyAugmented::AbstractGroup
       * for this method to be implemented.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyInverseTranspose(Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector& B,
                NOX::Abstract::MultiVector& X) const;

    private:

      //! Private to prohibit copying
      ComplexOperator(const ComplexOperator&);

      //! Private to prohibit copying
      ComplexOperator& operator = (const ComplexOperator&);

    protected:

      //! Pointer to group
      Teuchos::RCP<const LOCA::Hopf::MooreSpence::AbstractGroup> grpPtr;

      //! Complex frequency
      double omega;

    }; // Class ComplexOperator

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_COMPLEXOPERATOR_H
