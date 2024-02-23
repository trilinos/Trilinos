/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

#include "Teko_BlockInvDiagonalStrategy.hpp"
#include "Teko_InverseFactory.hpp"

namespace Teko {

InvFactoryDiagStrategy::InvFactoryDiagStrategy(const Teuchos::RCP<InverseFactory>& factory) {
  // only one factory to use!
  invDiagFact_.resize(1, factory);
  defaultInvFact_ = factory;
}

InvFactoryDiagStrategy::InvFactoryDiagStrategy(
    const std::vector<Teuchos::RCP<InverseFactory> >& factories,
    const Teuchos::RCP<InverseFactory>& defaultFact) {
  invDiagFact_ = factories;

  if (defaultFact == Teuchos::null)
    defaultInvFact_ = invDiagFact_[0];
  else
    defaultInvFact_ = defaultFact;
}

InvFactoryDiagStrategy::InvFactoryDiagStrategy(
    const std::vector<Teuchos::RCP<InverseFactory> >& inverseFactories,
    const std::vector<Teuchos::RCP<InverseFactory> >& preconditionerFactories,
    const Teuchos::RCP<InverseFactory>& defaultInverseFact,
    const Teuchos::RCP<InverseFactory>& defaultPreconditionerFact) {
  invDiagFact_  = inverseFactories;
  precDiagFact_ = preconditionerFactories;

  if (defaultInverseFact == Teuchos::null)
    defaultInvFact_ = invDiagFact_[0];
  else
    defaultInvFact_ = defaultInverseFact;
  defaultPrecFact_ = defaultPreconditionerFact;
}

/** returns an (approximate) inverse of the diagonal blocks of A
 * where A is closely related to the original source for invD0 and invD1
 */
void InvFactoryDiagStrategy::getInvD(const BlockedLinearOp& A, BlockPreconditionerState& state,
                                     std::vector<LinearOp>& invDiag) const {
  Teko_DEBUG_SCOPE("InvFactoryDiagStrategy::getInvD", 10);

  // loop over diagonals, build an inverse operator for each
  size_t diagCnt = A->productRange()->numBlocks();

  Teko_DEBUG_MSG("# diags = " << diagCnt << ", # inverses = " << invDiagFact_.size(), 6);

  const std::string opPrefix = "JacobiDiagOp";
  for (size_t i = 0; i < diagCnt; i++) {
    auto precFact = ((i < precDiagFact_.size()) && (!precDiagFact_[i].is_null()))
                        ? precDiagFact_[i]
                        : defaultPrecFact_;
    auto invFact  = (i < invDiagFact_.size()) ? invDiagFact_[i] : defaultInvFact_;
    invDiag.push_back(buildInverse(*invFact, precFact, getBlock(i, i, A), state, opPrefix, i));
  }
}

LinearOp InvFactoryDiagStrategy::buildInverse(const InverseFactory& invFact,
                                              Teuchos::RCP<InverseFactory>& precFact,
                                              const LinearOp& matrix,
                                              BlockPreconditionerState& state,
                                              const std::string& opPrefix, int i) const {
  std::stringstream ss;
  ss << opPrefix << "_" << i;

  ModifiableLinearOp& invOp  = state.getModifiableOp(ss.str());
  ModifiableLinearOp& precOp = state.getModifiableOp("prec_" + ss.str());

  if (precFact != Teuchos::null) {
    if (precOp == Teuchos::null) {
      precOp = precFact->buildInverse(matrix);
      state.addModifiableOp("prec_" + ss.str(), precOp);
    } else {
      Teko::rebuildInverse(*precFact, matrix, precOp);
    }
  }

  if (invOp == Teuchos::null)
    if (precOp.is_null())
      invOp = Teko::buildInverse(invFact, matrix);
    else
      invOp = Teko::buildInverse(invFact, matrix, precOp);
  else {
    if (precOp.is_null())
      Teko::rebuildInverse(invFact, matrix, invOp);
    else
      Teko::rebuildInverse(invFact, matrix, precOp, invOp);
  }

  return invOp;
}

}  // end namespace Teko
