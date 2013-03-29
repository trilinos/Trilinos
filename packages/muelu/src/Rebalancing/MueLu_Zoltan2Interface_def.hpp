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
#ifndef MUELU_ZOLTAN2INTERFACE_DEF_HPP
#define MUELU_ZOLTAN2INTERFACE_DEF_HPP

#include <sstream>

#include "MueLu_Zoltan2Interface_decl.hpp"
#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)

#include <Zoltan2_BasicCoordinateInput.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

#include <Teuchos_Utils.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

 template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
 RCP<const ParameterList> Zoltan2Interface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",           Teuchos::null, "Factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory of the coordinates");
    validParamList->set< RCP<const FactoryBase> >("number of partitions", Teuchos::null, "(advanced) Factory computing the number of partition.");

    return validParamList;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Zoltan2Interface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "Coordinates");
    Input(currentLevel, "number of partitions");

  } //DeclareInput()

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Zoltan2Interface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &level) const {
    FactoryMonitor m(*this, "Build", level);

    RCP<Matrix>         A = Get< RCP<Matrix> >(level, "A");
    GO           numParts = Get<GO>(level, "number of partitions");
    LocalOrdinal  blkSize = A->GetFixedBlockSize();

    if (!IsAvailable(level, "Coordinates")) {
      std::cout << GetFactory("Coordinates") << std::endl;

      level.print(*getOStream());
      throw Exceptions::RuntimeError("MueLu::ZoltanInterface::Build(): no coordinates available");
    }

    RCP<MultiVector> coords = Get< RCP<MultiVector> >(level, "Coordinates");
    size_t              dim = coords->getNumVectors();

    RCP<const Map> map = coords->getMap();

    if (numParts == 1) {
      // Single processor, decomposition is trivial: all zeros
      RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(map, true);
      Set(level, "Partition", decomposition);
      return;
    }

    GO numElements = map->getNodeNumElements();
    std::vector<const SC*> values(dim), weights(1);
    std::vector<int>       strides;

    for (size_t k = 0; k < dim; k++)
      values[k] = coords->getData(k).get();

    Array<SC> numEntriesPerRow(numElements);
    for (LO i = 0; i < numElements; i++)
      for (LO j = 0; j < blkSize; j++)
        numEntriesPerRow[i] += A->getNumEntriesInLocalRow(i*blkSize+j);
    weights[0] = numEntriesPerRow.getRawPtr();

    Teuchos::ParameterList params;

    bool usePQJagged = true;
    if (usePQJagged) {
      params.set("algorithm",             "multijagged");

      // Find representation of numParts as a product of dim factors.
      // This is needed for pqJagged algorithm, as it splits each dimension independently
      //
      // We would like the factors to be as close to each other as possible
      // While there is some linear search, the worst case for 1billion elements
      // would produce two linear searches of length sqrt(1 billion) = 10000, which
      // should not be absolutely terrible
      // For now, we ignore the fact that we might get smth like (1,1,937)
      std::ostringstream pqParts;
      if (dim == 1)
        pqParts << numParts;

      else if (dim == 2) {
        GO p = Teuchos::as<GO>(sqrt(numParts));
        while (numParts % p)
          p--;
        pqParts << p << "," << numParts/p;

      } else if (dim == 3) {
        GO p = Teuchos::as<GO>(pow(numParts,1./3)), q;
        while (numParts % p)
          p--;
        q = Teuchos::as<GO>(sqrt(numParts / p));
        while (numParts % (p*q))
          q--;
        pqParts << p << "," << q << "," << numParts/(p*q);
      }

      params.set("pqParts",               pqParts.str());

    } else {
      params.set("algorithm",             "rcb");
      params.set("num_global_parts",      numParts);
    }

    params.set("partitioning_approach",   "partition");
    // params.set("compute_metrics",         "true");

    typedef Zoltan2::BasicCoordinateInput<Zoltan2::BasicUserTypes<SC,GO,LO,GO> > InputAdapterType;
    typedef Zoltan2::PartitioningProblem<InputAdapterType> ProblemType;

    InputAdapterType adapter(numElements, map->getNodeElementList().getRawPtr(), values, strides, weights, strides);

    const Teuchos::MpiComm<int>& comm = static_cast<const Teuchos::MpiComm<int>& >(*map->getComm());
    RCP<ProblemType> problem(new ProblemType(&adapter, &params, *comm.getRawMpiComm()));

    problem->solve();

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(A->getRowMap(), false);
    ArrayRCP<GO>                      decompEntries = decomposition->getDataNonConst(0);

    const zoltan2_partId_t * parts = problem->getSolution().getPartList();

    for (GO localID = 0; localID < numElements; localID++) {
      int partNum = parts[localID];

      for (LO j = 0; j < blkSize; j++)
        decompEntries[localID*blkSize + j] = partNum;
    }

    Set(level, "Partition", decomposition);

  } //Build()

} //namespace MueLu

#endif //if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)

#endif // MUELU_ZOLTAN2INTERFACE_DEF_HPP
