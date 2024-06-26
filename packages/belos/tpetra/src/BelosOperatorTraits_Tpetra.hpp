// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSOPERATORTRAITS_TPETRA_HPP
#define BELOSOPERATORTRAITS_TPETRA_HPP

/// \file BelosMultiVecTraits_Tpetra.hpp
/// \brief Partial specialization of Belos::MultiVecTraits for Tpetra objects.

#include "BelosOperatorTraits.hpp"
#include "Tpetra_Operator.hpp"

namespace Belos {

  //! Partial specialization of OperatorTraits for Tpetra objects.
  template <class Scalar, class LO, class GO, class Node>
  class OperatorTraits<Scalar,
                       ::Tpetra::MultiVector<Scalar,LO,GO,Node>,
                       ::Tpetra::Operator<Scalar,LO,GO,Node> >
  {
  public:
    static void
    Apply (const ::Tpetra::Operator<Scalar,LO,GO,Node>& Op,
           const ::Tpetra::MultiVector<Scalar,LO,GO,Node>& X,
           ::Tpetra::MultiVector<Scalar,LO,GO,Node>& Y,
           const ETrans trans = NOTRANS)
    {
      Teuchos::ETransp teuchosTrans = Teuchos::NO_TRANS;
      if (trans == NOTRANS) {
        teuchosTrans = Teuchos::NO_TRANS;
      } else if (trans == TRANS) {
        teuchosTrans = Teuchos::TRANS;
      } else if (trans == CONJTRANS) {
        teuchosTrans = Teuchos::CONJ_TRANS;
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::invalid_argument, "Belos::OperatorTraits::Apply: Invalid "
          "'trans' value " << trans << ".  Valid values are NOTRANS=" << NOTRANS
          << ", TRANS=" << TRANS << ", and CONJTRANS=" << CONJTRANS << ".");
      }
      Op.apply (X, Y, teuchosTrans);
    }

    static bool
    HasApplyTranspose (const ::Tpetra::Operator<Scalar,LO,GO,Node>& Op)
    {
      return Op.hasTransposeApply ();
    }
  };

} // namespace Belos

#endif // BELOSOPERATORTRAITS_TPETRA_HPP
