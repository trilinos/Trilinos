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
#ifndef PACKAGES_MUELU_SRC_REBALANCING_MUELU_CLONEREPARTITIONINTERFACE_DEF_HPP_
#define PACKAGES_MUELU_SRC_REBALANCING_MUELU_CLONEREPARTITIONINTERFACE_DEF_HPP_

#include "MueLu_CloneRepartitionInterface_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"


namespace MueLu {

 template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
 Teuchos::RCP<const ParameterList> CloneRepartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    Teuchos::RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< Teuchos::RCP<const FactoryBase> >("A",         Teuchos::null, "Factory of the matrix A");
    validParamList->set< Teuchos::RCP<const FactoryBase> >("Partition", Teuchos::null, "Factory generating the Partition array (e.g. ZoltanInterface)");

    return validParamList;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CloneRepartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level & currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "Partition");
  } //DeclareInput()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CloneRepartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);
    currentLevel.print(GetOStream(Statistics0,0));
    // extract blocked operator A from current level
    Teuchos::RCP<Matrix> A = Get< Teuchos::RCP<Matrix> >     (currentLevel, "A");
    Teuchos::RCP<const Teuchos::Comm< int > > comm = A->getRowMap()->getComm();
    //const int myRank = comm->getRank();

    // number of Partitions only used for a shortcut.
    GO numPartitions = 0;
    if (currentLevel.IsAvailable("number of partitions")) {
      numPartitions = currentLevel.Get<GO>("number of partitions");
      GetOStream(Warnings0) << "Using user-provided \"number of partitions\", the performance is unknown" << std::endl;

    }

    // ======================================================================================================
    // Construct decomposition vector
    // ======================================================================================================
    RCP<GOVector> decomposition = Teuchos::null;
    if (numPartitions == 1) {
      // Trivial case: decomposition is the trivial one, all zeros. We skip the call to Zoltan_Interface
      // (this is mostly done to avoid extra output messages, as even if we didn't skip there is a shortcut
      // in Zoltan[12]Interface).
      GetOStream(Warnings0) << "Only one partition: Skip call to the repartitioner." << std::endl;
      decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), true);
      Set(currentLevel, "Partition", decomposition);
      return;
    }

    // extract decomposition vector
    decomposition = Get<RCP<GOVector> >(currentLevel, "Partition");
    ArrayRCP<const GO> decompEntries = decomposition->getData(0);

    if (decomposition.is_null()) {
      GetOStream(Warnings0) << "No repartitioning necessary: partitions were left unchanged by the repartitioner" << std::endl;
      Set<RCP<const Import> >(currentLevel, "Importer", Teuchos::null);
      return;
    }



    /*std::vector<int> amActive  = std::vector<int>(comm->getSize(),0);
    std::vector<int> areActive = std::vector<int>(comm->getSize(),0);
    if() amActive[myRank] = 1;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX,comm->getSize(),&amActive[0],&areActive[0]);*/



    // plausibility check!

    size_t inLocalLength  = decomposition->getLocalLength();
    size_t outLocalLength = A->getRowMap()->getNodeNumElements();
    Scalar ratioLocal = Teuchos::as<Scalar>(inLocalLength) / Teuchos::as<Scalar>(outLocalLength);
    Xpetra::global_size_t inGlobalLength  = decomposition->getGlobalLength();
    Xpetra::global_size_t outGlobalLength = A->getRowMap()->getGlobalNumElements();
    Scalar ratioGlobal = Teuchos::as<Scalar>(inGlobalLength) / Teuchos::as<Scalar>(outGlobalLength);

    TEUCHOS_TEST_FOR_EXCEPTION(ratioLocal != ratioGlobal, MueLu::Exceptions::RuntimeError,"CloneRepartitionInterface: ratio of input and output vector lengths is not consistent with global information.");
    TEUCHOS_TEST_FOR_EXCEPTION(inLocalLength % outLocalLength != 0, MueLu::Exceptions::RuntimeError,"CloneRepartitionInterface: ratio of input and output vector lengths is not consistent with global information.");
    TEUCHOS_TEST_FOR_EXCEPTION(inGlobalLength % outGlobalLength != 0, MueLu::Exceptions::RuntimeError,"CloneRepartitionInterface: ratio of input and output vector lengths is not consistent with global information.");

    // create new decomposition vector
    Teuchos::RCP<GOVector> ret = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), false);
    ArrayRCP<GO> retDecompEntries = ret->getDataNonConst(0);

    // fill decomposition vector
    LO ratio = Teuchos::as<LO>(inLocalLength / outLocalLength);
    for(LO i = 0; i<Teuchos::as<LO>(ret->getMap()->getNodeNumElements()); i++) {
      retDecompEntries[i] = Teuchos::as<GO>(decompEntries[i*ratio]);
    }
    Set(currentLevel, "Partition", ret);
  } //Build()



} //namespace MueLu

#endif /* PACKAGES_MUELU_SRC_REBALANCING_MUELU_CLONEREPARTITIONINTERFACE_DEF_HPP_ */
