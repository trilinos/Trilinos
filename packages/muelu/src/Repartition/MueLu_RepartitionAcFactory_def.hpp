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
#ifndef MUELU_REPARTITIONACFACTORY_DEF_HPP
#define MUELU_REPARTITIONACFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_RepartitionAcFactory_decl.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(coarseLevel, "P");
    Input(coarseLevel, "Importer");
    coarseLevel.DeclareInput("A", GetFactory("P").get(), this); //FIXME hack
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Computing Ac", coarseLevel);

    RCP<Matrix> P = Get< RCP<Matrix> >(coarseLevel, "P");

    RCP<Matrix> Ac;

    if ( coarseLevel.IsAvailable("A",GetFactory("P").get()) && IsAvailable(coarseLevel, "Importer") ) {

      SubFactoryMonitor subM(*this, "Rebalancing existing Ac", coarseLevel);

      Ac = coarseLevel.Get< RCP<Matrix> >("A", GetFactory("P").get());
      RCP<Matrix> newAc = MatrixFactory::Build(P->getDomainMap(), Ac->getGlobalMaxNumRowEntries());
      RCP<CrsMatrixWrap> crsOp = rcp_dynamic_cast<CrsMatrixWrap>(newAc);
      RCP<CrsMatrix> crsMtx = crsOp->getCrsMatrix();
      RCP<CrsMatrixWrap> origOp = rcp_dynamic_cast<CrsMatrixWrap>(Ac);
      RCP<CrsMatrix> origMtx = origOp->getCrsMatrix();
      RCP<const Import> permImporter = Get< RCP<const Import> >(coarseLevel, "Importer");
      crsMtx->doImport(*origMtx, *permImporter,Xpetra::INSERT);
      crsMtx = Teuchos::null;

      //TODO add plausibility check

      newAc->fillComplete(P->getDomainMap(), P->getDomainMap());
      Ac = newAc;  // TODO: rename variables to avoid this swap operation

      GetOStream(Statistics0, 0) << RAPFactory::PrintMatrixInfo(*Ac, "Ac (rebalanced)");
      GetOStream(Statistics0, 0) << RAPFactory::PrintLoadBalancingInfo(*Ac, "Ac (rebalanced)");

    } else if (coarseLevel.IsAvailable("A",GetFactory("P").get())) {

      // Ac already built by the load balancing process and no load balancing needed
      SubFactoryMonitor subM(*this, "Ac already computed", coarseLevel);
      Ac = coarseLevel.Get< RCP<Matrix> >("A", GetFactory("P").get());

    } else {

      TEUCHOS_TEST_FOR_EXCEPT(true);

    }

    Set(coarseLevel, "A", Ac);

  } //Build()

} //namespace MueLu

#define MUELU_REPARTITIONACFACTORY_SHORT
#endif // MUELU_REPARTITIONACFACTORY_DEF_HPP
