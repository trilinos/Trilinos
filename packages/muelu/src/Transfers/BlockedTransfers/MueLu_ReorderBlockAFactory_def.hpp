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

#include <Xpetra_BlockedCrsMatrix.hpp>
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
    std::string reorderStr = Teuchos::as<std::string>(pL.get<std::string>("Reorder Type"));

    RCP<Matrix>           Ain = Get<RCP<Matrix> >(currentLevel, "A");
    RCP<BlockedCrsMatrix> A   = rcp_dynamic_cast<BlockedCrsMatrix>(Ain);

    TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(),     Exceptions::BadCast,      "Input matrix A is not a BlockedCrsMatrix.");

    //Teuchos::RCP<Block::BlockReorderManager> brm = blockedReorderFromString("[[0 1] 2]");


    //RCP<Matrix> Op = A->getMatrix(row, col);

    // build new map extractors containing sub blocks
    //RCP<const MapExtractor> rangeMapExtractor  = A->getRangeMapExtractor();
    //RCP<const MapExtractor> domainMapExtractor = A->getDomainMapExtractor();


    /*GetOStream(Statistics1) << "A(" << row << "," << col << ") has strided maps:"
        << "\n  range  map fixed block size = " << srangeMap ->getFixedBlockSize() << ", strided block id = " << srangeMap ->getStridedBlockId()
        << "\n  domain map fixed block size = " << sdomainMap->getFixedBlockSize() << ", strided block id = " << sdomainMap->getStridedBlockId() << std::endl;

    // TODO do we really need that? we moved the code to getMatrix...
    if (Op->IsView("stridedMaps") == true)
      Op->RemoveView("stridedMaps");
    Op->CreateView("stridedMaps", srangeMap, sdomainMap);

    TEUCHOS_TEST_FOR_EXCEPTION(Op->IsView("stridedMaps") == false, Exceptions::RuntimeError, "Failed to set \"stridedMaps\" view.");
    */
    //currentLevel.Set("A", Op, this);
  }

} // namespace MueLu

#endif /* MUELU_REORDERBLOCKAFACTORY_DEF_HPP_ */

