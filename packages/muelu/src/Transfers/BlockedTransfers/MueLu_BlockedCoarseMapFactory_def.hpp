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
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_BlockedCoarseMapFactory_def.hpp
 *
 *  Created on: Oct 16, 2012
 *      Author: tobias
 */

#ifndef MUELU_BLOCKEDCOARSEMAPFACTORY_DEF_HPP_
#define MUELU_BLOCKEDCOARSEMAPFACTORY_DEF_HPP_

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_BlockedCoarseMapFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_Factory.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedCoarseMapFactory()
  {  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~BlockedCoarseMapFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Generating factory for aggregates.");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",  Teuchos::null, "Generating factory for null space.");
    validParamList->set< RCP<const FactoryBase> >("CoarseMap",  Teuchos::null, "Generating factory of previous coarse map. (must be set by user!).");

    // do we need this?
    validParamList->set< std::string  >("Striding info", "{}", "Striding information");
    validParamList->set< LocalOrdinal >("Strided block id", -1, "Strided block id");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    this->Input(currentLevel, "Aggregates");
    this->Input(currentLevel, "Nullspace");

    // Get CoarseMap from previously defined block
    RCP<const FactoryBase> prevCoarseMapFact = this->GetFactory("CoarseMap");
    TEUCHOS_TEST_FOR_EXCEPTION(prevCoarseMapFact==Teuchos::null, Exceptions::RuntimeError, "MueLu::BlockedCoarseMapFactory::getDomainMapOffset: user did not specify CoarseMap of previous block. Do not forget to set the CoarseMap factory.");
    currentLevel.DeclareInput("CoarseMap", prevCoarseMapFact.get(), this); // --
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "BlockedCoarseMap factory", currentLevel);

    RCP<const FactoryBase> prevCoarseMapFact = this->GetFactory("CoarseMap");
    RCP<const Map> subPDomainMap = currentLevel.Get<RCP<const Map> >("CoarseMap", prevCoarseMapFact.get() /*prevCoarseMapFact_.get()*/);

    GlobalOrdinal maxGlobalIndex = subPDomainMap->getMaxAllGlobalIndex();

    RCP<Aggregates> aggregates = Factory::Get< RCP<Aggregates> >(currentLevel, "Aggregates");
    GlobalOrdinal numAggs = aggregates->GetNumAggregates();

    // extract communicator
    RCP<const Teuchos::Comm<int> > comm = aggregates->GetMap()->getComm();

    // determine nullspace dimension
    RCP<MultiVector> nullspace  = Factory::Get< RCP<MultiVector> >(currentLevel, "Nullspace");
    const size_t NSDim = nullspace->getNumVectors();

    LocalOrdinal stridedBlockId = CoarseMapFactory::getStridedBlockId();

    // check for consistency of striding information with NSDim and nCoarseDofs
    if( stridedBlockId== -1 ) {
      // this means we have no real strided map but only a block map with constant blockSize "NSDim"
      TEUCHOS_TEST_FOR_EXCEPTION(CoarseMapFactory::stridingInfo_.size() > 1, Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): stridingInfo_.size() but must be one");
      CoarseMapFactory::stridingInfo_.clear();
      CoarseMapFactory::stridingInfo_.push_back(NSDim);
      TEUCHOS_TEST_FOR_EXCEPTION(CoarseMapFactory::stridingInfo_.size() != 1, Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): stridingInfo_.size() but must be one");
    } else {
      // stridedBlockId_ > -1, set by user
      TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockId > Teuchos::as<LO>(CoarseMapFactory::stridingInfo_.size() - 1) , Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): it is stridingInfo_.size() <= stridedBlockId_. error.");
      size_t stridedBlockSize = CoarseMapFactory::stridingInfo_[stridedBlockId];
      TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockSize != NSDim , Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): dimension of strided block != NSDim. error.");
    }

    CoarseMapFactory::GetOStream(Statistics2) << "domainGIDOffset: " << maxGlobalIndex + 1 << " block size: " << CoarseMapFactory::getFixedBlockSize() << " stridedBlockId: " << stridedBlockId << std::endl;

    // number of coarse level dofs (fixed by number of aggregates and blocksize data)
    GlobalOrdinal nCoarseDofs = numAggs * CoarseMapFactory::getFixedBlockSize();
    GlobalOrdinal indexBase   = aggregates->GetMap()->getIndexBase();

    RCP<const Map> coarseMap = StridedMapFactory::Build(aggregates->GetMap()->lib(),
        Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        nCoarseDofs,
        indexBase,
        CoarseMapFactory::stridingInfo_,
        comm,
        stridedBlockId,
        maxGlobalIndex + 1);

    this->Set(currentLevel, "CoarseMap", coarseMap);
  } // Build

} //namespace MueLu

#endif /* MUELU_BLOCKEDCOARSEMAPFACTORY_DEF_HPP_ */
