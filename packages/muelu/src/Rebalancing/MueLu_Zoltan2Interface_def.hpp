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
#ifndef MUELU_ZOLTAN2INTERFACE_DEF_HPP
#define MUELU_ZOLTAN2INTERFACE_DEF_HPP

#include <sstream>
#include <set>

#include "MueLu_Zoltan2Interface_decl.hpp"
#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)

#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

#include <Teuchos_Utils.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

 template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
 Zoltan2Interface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Zoltan2Interface() {
    defaultZoltan2Params = rcp(new ParameterList());
    defaultZoltan2Params->set("algorithm",             "multijagged");
    defaultZoltan2Params->set("partitioning_approach", "partition");
 }

 template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
 RCP<const ParameterList> Zoltan2Interface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >   ("A",                      Teuchos::null, "Factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >   ("Coordinates",            Teuchos::null, "Factory of the coordinates");
    validParamList->set< int >                      ("rowWeight",                          0, "Default weight to rows (total weight = nnz + rowWeight");
    validParamList->set< RCP<const ParameterList> > ("ParameterList",          Teuchos::null, "Zoltan2 parameters");

    return validParamList;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Zoltan2Interface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "Coordinates");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Zoltan2Interface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& level) const {
    FactoryMonitor m(*this, "Build", level);

    RCP<Matrix>      A        = Get< RCP<Matrix> >     (level, "A");
    RCP<const Map>   rowMap   = A->getRowMap();

    RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > coords   = Get< RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > >(level, "Coordinates");
    RCP<const Map>   map      = coords->getMap();

    GO               numParts = level.Get<GO>("number of partitions");

    size_t dim       = coords->getNumVectors();
    LO     blkSize   = A->GetFixedBlockSize();

    // Check that the number of local coordinates is consistent with the #rows in A
    TEUCHOS_TEST_FOR_EXCEPTION(rowMap->getNodeNumElements()/blkSize != coords->getLocalLength(), Exceptions::Incompatible,
                               "Coordinate vector length (" + toString(coords->getLocalLength()) << " is incompatible with number of block rows in A ("
                               + toString(rowMap->getNodeNumElements()/blkSize) + "The vector length should be the same as the number of mesh points.");
#ifdef HAVE_MUELU_DEBUG
    GO indexBase = rowMap->getIndexBase();
    GetOStream(Runtime0) << "Checking consistence of row and coordinates maps" << std::endl;
    // Make sure that logical blocks in row map coincide with logical nodes in coordinates map
    ArrayView<const GO> rowElements    = rowMap->getNodeElementList();
    ArrayView<const GO> coordsElements = map   ->getNodeElementList();
    for (LO i = 0; i < Teuchos::as<LO>(map->getNodeNumElements()); i++)
      TEUCHOS_TEST_FOR_EXCEPTION((coordsElements[i]-indexBase)*blkSize + indexBase != rowElements[i*blkSize],
                                 Exceptions::RuntimeError, "i = " << i << ", coords GID = " << coordsElements[i]
                                 << ", row GID = " << rowElements[i*blkSize] << ", blkSize = " << blkSize << std::endl);
#endif

    if (numParts == 1) {
      // Single processor, decomposition is trivial: all zeros
      RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
      Set(level, "Partition", decomposition);
      return;
    }

    GO numElements = map->getNodeNumElements();
    std::vector<const double*> values(dim), weights(1);
    std::vector<int>       strides;

    for (size_t k = 0; k < dim; k++)
      values[k] = coords->getData(k).get();

    const ParameterList& pL = GetParameterList();
    int rowWeight = pL.get<int>("rowWeight");
    GetOStream(Runtime0) << "Using weights formula: nnz + " << rowWeight << std::endl;

    Array<double> weightsPerRow(numElements);
    for (LO i = 0; i < numElements; i++) {
      weightsPerRow[i] = Teuchos::ScalarTraits<SC>::zero();
      for (LO j = 0; j < blkSize; j++) {
        weightsPerRow[i] += A->getNumEntriesInLocalRow(i*blkSize+j);
        // Zoltan2 pqJagged gets as good partitioning as Zoltan RCB in terms of nnz
        // but Zoltan also gets a good partioning in rows, which sometimes does not
        // happen for Zoltan2. So here is an attempt to get a better row partitioning
        // without significantly screwing up nnz partitioning
        // NOTE: no good heuristic here, the value was chosen almost randomly
        weightsPerRow[i] += rowWeight;
      }
    }
    weights[0] = weightsPerRow.getRawPtr();

    RCP<const ParameterList> providedList = pL.get<RCP<const ParameterList> >("ParameterList");
    ParameterList Zoltan2Params;
    if (providedList != Teuchos::null)
      Zoltan2Params = *providedList;

    // Merge defalt Zoltan2 parameters with user provided
    // If default and user parameters contain the same parameter name, user one is always preferred
    for (ParameterList::ConstIterator param = defaultZoltan2Params->begin(); param != defaultZoltan2Params->end(); param++) {
      const std::string& pName = defaultZoltan2Params->name(param);
      if (!Zoltan2Params.isParameter(pName))
        Zoltan2Params.set(pName, defaultZoltan2Params->get<std::string>(pName));
    }
    Zoltan2Params.set("num_global_parts", Teuchos::as<int>(numParts));

    GetOStream(Runtime0) << "Zoltan2 parameters:" << std::endl << "----------" << std::endl << Zoltan2Params << "----------" << std::endl;

    const std::string& algo = Zoltan2Params.get<std::string>("algorithm");
    TEUCHOS_TEST_FOR_EXCEPTION(algo != "multijagged" &&
                               algo != "rcb",
                               Exceptions::RuntimeError, "Unknown partitioning algorithm: \"" << algo << "\"");

    typedef Zoltan2::BasicVectorAdapter<Zoltan2::BasicUserTypes<double,GO,LO,GO> > InputAdapterType;
    typedef Zoltan2::PartitioningProblem<InputAdapterType> ProblemType;

    InputAdapterType adapter(numElements, map->getNodeElementList().getRawPtr(), values, strides, weights, strides);

    RCP<const Teuchos::MpiComm<int> >            dupMpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(rowMap->getComm()->duplicate());
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > zoltanComm = dupMpiComm->getRawMpiComm();

    RCP<ProblemType> problem(new ProblemType(&adapter, &Zoltan2Params, (*zoltanComm)()));

    {
      SubFactoryMonitor m1(*this, "Zoltan2 " + toString(algo), level);
      problem->solve();
    }

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(rowMap, false);
    ArrayRCP<GO>                      decompEntries = decomposition->getDataNonConst(0);

    const typename InputAdapterType::part_t * parts = problem->getSolution().getPartListView();

    for (GO i = 0; i < numElements; i++) {
      int partNum = parts[i];

      for (LO j = 0; j < blkSize; j++)
        decompEntries[i*blkSize + j] = partNum;
    }

    Set(level, "Partition", decomposition);

  } //Build()

} //namespace MueLu

#endif //if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)

#endif // MUELU_ZOLTAN2INTERFACE_DEF_HPP
