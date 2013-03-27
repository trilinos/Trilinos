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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_REBALANCEACFACTORY_DEF_HPP
#define MUELU_REBALANCEACFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_RebalanceAcFactory_decl.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("A",         Teuchos::null, "Generating factory of the matrix A for rebalancing");
    validParamList->set< RCP<const FactoryBase> >("Importer",  Teuchos::null, "Generating factory of the importer");
    // The value of "useSubcomm" paramter here must be the same as in RebalanceTransferFactory
    validParamList->set< bool >                  ("useSubcomm",         true, "Construct subcommunicators");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(coarseLevel, "A");
    Input(coarseLevel, "Importer");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Computing Ac", coarseLevel);

    RCP<Matrix> originalAc = Get< RCP<Matrix> >(coarseLevel, "A");

    RCP<const Import> rebalanceImporter = Get< RCP<const Import> >(coarseLevel, "Importer");

    if (rebalanceImporter != Teuchos::null) {

      RCP<Matrix> rebalancedAc;
      {
        SubFactoryMonitor subM(*this, "Rebalancing existing Ac", coarseLevel);

        RCP<const Map> targetMap = rebalanceImporter->getTargetMap();	
        rebalancedAc = MatrixFactory::Build(targetMap, originalAc->getGlobalMaxNumRowEntries());

        rebalancedAc->doImport(*originalAc, *rebalanceImporter, Xpetra::INSERT);

        const ParameterList & pL = GetParameterList();
        if (pL.get<bool>("useSubcomm") == true) {
          // replace full communicator with a subcommunicator
          GetOStream(Runtime0,0) << "Replacing maps with a subcommunicator" << std::endl;

          RCP<const Map> reducedMap = targetMap->removeEmptyProcesses();
          rebalancedAc->removeEmptyProcessesInPlace(reducedMap);
          if (!reducedMap.is_null()) {
            // we own some part of the new matrix
            rebalancedAc->fillComplete();
            rebalancedAc->SetFixedBlockSize(originalAc->GetFixedBlockSize());
          } else {
            rebalancedAc = Teuchos::null;
          }

        } else {
          rebalancedAc->fillComplete();
          rebalancedAc->SetFixedBlockSize(originalAc->GetFixedBlockSize());
        }

#ifdef OLD_AND_BUSTED
	rebalancedAc = MatrixFactory::Build(originalAc,*rebalanceImporter,targetMap,targetMap);
        RCP<const Map> targetMap = rebalanceImporter->getTargetMap();
        rebalancedAc->SetFixedBlockSize(originalAc->GetFixedBlockSize());
#endif
        Set(coarseLevel, "A", rebalancedAc);
      }

      if (!rebalancedAc.is_null()) {
        GetOStream(Statistics0, 0) << RAPFactory::PrintMatrixInfo       (*rebalancedAc, "Ac (rebalanced)");
        GetOStream(Statistics0, 0) << RAPFactory::PrintLoadBalancingInfo(*rebalancedAc, "Ac (rebalanced)");
      }

    } else {

      // Ac already built by the load balancing process and no load balancing needed
      GetOStream(Warnings0, 0) << "No rebalancing" << std::endl;
      GetOStream(Warnings0, 0) << "Jamming A into Level " << coarseLevel.GetLevelID() << " w/ generating factory "
                               << this << std::endl;

      Set(coarseLevel, "A", originalAc);
    }

  } //Build()

} //namespace MueLu

#endif // MUELU_REBALANCEACFACTORY_DEF_HPP
