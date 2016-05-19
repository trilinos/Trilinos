// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef MUELU_REORDERBLOCKAFACTORY_DEF_HPP_
#define MUELU_REORDERBLOCKAFACTORY_DEF_HPP_


#include "MueLu_ReorderBlockAFactory_decl.hpp"

#include <Xpetra_BlockReorderManager.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> ReorderBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",          MueLu::NoFactory::getRCP(), "Generating factory for A.");

    validParamList->set< std::string  >          ("Reorder Type",    "", "String describing the reordering of blocks");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ReorderBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ReorderBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    const ParameterList& pL = GetParameterList();
    std::string reorderStr = pL.get<std::string>("Reorder Type");

    RCP<Matrix>           Ain = Get<RCP<Matrix> >(currentLevel, "A");
    RCP<BlockedCrsMatrix> A   = rcp_dynamic_cast<BlockedCrsMatrix>(Ain);

    TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(),     Exceptions::BadCast,      "Input matrix A is not a BlockedCrsMatrix.");

    Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString(reorderStr);

    Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop =
        Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(Xpetra::buildReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(brm, A));

    TEUCHOS_TEST_FOR_EXCEPTION(brop.is_null(),     Exceptions::RuntimeError,      "Block reordering of " << A->Rows() << "x" << A->Cols() << " blocked operator failed.");

    GetOStream(Statistics1) << "Reordering A using " << brm->toString() << " block gives a " << brop->Rows() << "x" << brop->Cols() << " blocked operators" << std::endl;
    GetOStream(Debug) << "Reordered operator has " << brop->getRangeMap()->getGlobalNumElements() << " rows and " << brop->getDomainMap()->getGlobalNumElements() << " columns" << std::endl;
    GetOStream(Debug) << "Reordered operator: Use of Thyra style gids = " << brop->getRangeMapExtractor()->getThyraMode() << std::endl;

    // get rid of const (we expect non-const operators stored in Level)
    Teuchos::RCP<ReorderedBlockedCrsMatrix> bret =
        Teuchos::rcp_const_cast<ReorderedBlockedCrsMatrix>(brop);

    // TODO strided maps
    // blocked operators do not have strided maps (information could be misleading?)
    //Op->CreateView("stridedMaps", srangeMap, sdomainMap);

    currentLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bret), this);
  }

} // namespace MueLu

#endif /* MUELU_REORDERBLOCKAFACTORY_DEF_HPP_ */

