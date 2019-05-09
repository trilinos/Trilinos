//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

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
