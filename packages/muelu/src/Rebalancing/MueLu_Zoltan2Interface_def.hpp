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

#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

#include <Teuchos_Utils.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Zoltan2Interface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Zoltan2Interface() {
  defaultZoltan2Params = rcp(new ParameterList());
  defaultZoltan2Params->set("algorithm", "multijagged");
  defaultZoltan2Params->set("partitioning_approach", "partition");

  // Improve scaling for communication bound algorithms by premigrating
  // coordinates to a subset of processors.
  // For more information, see Github issue #1538
  defaultZoltan2Params->set("mj_premigration_option", 1);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> Zoltan2Interface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("number of partitions", Teuchos::null, "Instance of RepartitionHeuristicFactory.");
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory of the coordinates");
  validParamList->set<RCP<const ParameterList> >("ParameterList", Teuchos::null, "Zoltan2 parameters");
  validParamList->set<RCP<const FactoryBase> >("repartition: heuristic target rows per process", Teuchos::null, "Factory for number of rows per process to use with MultiJagged");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Zoltan2Interface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "number of partitions");
  const ParameterList& pL = GetParameterList();
  // We do this dance, because we don't want "ParameterList" to be marked as used.
  // Is there a better way?
  Teuchos::ParameterEntry entry                  = pL.getEntry("ParameterList");
  RCP<const Teuchos::ParameterList> providedList = Teuchos::any_cast<RCP<const Teuchos::ParameterList> >(entry.getAny(false));
  if (providedList != Teuchos::null && providedList->isType<std::string>("algorithm")) {
    const std::string algo = providedList->get<std::string>("algorithm");
    if (algo == "multijagged") {
      Input(currentLevel, "Coordinates");
      Input(currentLevel, "repartition: heuristic target rows per process");
    } else if (algo == "rcb") {
      Input(currentLevel, "Coordinates");
    }
  } else {
    Input(currentLevel, "repartition: heuristic target rows per process");
    Input(currentLevel, "Coordinates");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Zoltan2Interface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& level) const {
  FactoryMonitor m(*this, "Build", level);

  typedef typename Teuchos::ScalarTraits<SC>::coordinateType real_type;
  typedef typename Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  RCP<Matrix> A         = Get<RCP<Matrix> >(level, "A");
  RCP<const Map> rowMap = A->getRowMap();
  LO blkSize            = A->GetFixedBlockSize();

  int numParts = Get<int>(level, "number of partitions");
  if (numParts == 1 || numParts == -1) {
    // Single processor, decomposition is trivial: all zeros
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
    Set(level, "Partition", decomposition);
    return;
  } /* else if (numParts == -1) {
     // No repartitioning
     RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Teuchos::null; //Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
     //decomposition->putScalar(Teuchos::as<Scalar>(rowMap->getComm()->getRank()));
     Set(level, "Partition", decomposition);
     return;
   }*/

  const ParameterList& pL = GetParameterList();

  RCP<const ParameterList> providedList = pL.get<RCP<const ParameterList> >("ParameterList");
  ParameterList Zoltan2Params;
  if (providedList != Teuchos::null)
    Zoltan2Params = *providedList;

  // Merge defalt Zoltan2 parameters with user provided
  // If default and user parameters contain the same parameter name, user one is always preferred
  for (ParameterList::ConstIterator param = defaultZoltan2Params->begin(); param != defaultZoltan2Params->end(); param++) {
    const std::string& pName = defaultZoltan2Params->name(param);
    if (!Zoltan2Params.isParameter(pName))
      Zoltan2Params.setEntry(pName, defaultZoltan2Params->getEntry(pName));
  }
  Zoltan2Params.set("num_global_parts", Teuchos::as<int>(numParts));

  GetOStream(Runtime0) << "Zoltan2 parameters:\n----------\n"
                       << Zoltan2Params << "----------" << std::endl;

  const std::string& algo = Zoltan2Params.get<std::string>("algorithm");

  if (algo == "multijagged" || algo == "rcb") {
    RCP<RealValuedMultiVector> coords = Get<RCP<RealValuedMultiVector> >(level, "Coordinates");
    RCP<const Map> map                = coords->getMap();
    GO numElements                    = map->getLocalNumElements();

    // Check that the number of local coordinates is consistent with the #rows in A
    TEUCHOS_TEST_FOR_EXCEPTION(rowMap->getLocalNumElements() / blkSize != coords->getLocalLength(), Exceptions::Incompatible,
                               "Coordinate vector length (" + toString(coords->getLocalLength()) << " is incompatible with number of block rows in A (" + toString(rowMap->getLocalNumElements() / blkSize) + "The vector length should be the same as the number of mesh points.");
#ifdef HAVE_MUELU_DEBUG
    GO indexBase = rowMap->getIndexBase();
    GetOStream(Runtime0) << "Checking consistence of row and coordinates maps" << std::endl;
    // Make sure that logical blocks in row map coincide with logical nodes in coordinates map
    ArrayView<const GO> rowElements    = rowMap->getLocalElementList();
    ArrayView<const GO> coordsElements = map->getLocalElementList();
    for (LO i = 0; i < Teuchos::as<LO>(numElements); i++)
      TEUCHOS_TEST_FOR_EXCEPTION((coordsElements[i] - indexBase) * blkSize + indexBase != rowElements[i * blkSize],
                                 Exceptions::RuntimeError, "i = " << i << ", coords GID = " << coordsElements[i] << ", row GID = " << rowElements[i * blkSize] << ", blkSize = " << blkSize << std::endl);
#endif

    typedef Zoltan2::XpetraMultiVectorAdapter<RealValuedMultiVector> InputAdapterType;
    typedef Zoltan2::PartitioningProblem<InputAdapterType> ProblemType;

    Array<real_type> weightsPerRow(numElements);
    for (LO i = 0; i < numElements; i++) {
      weightsPerRow[i] = 0.0;

      for (LO j = 0; j < blkSize; j++) {
        weightsPerRow[i] += A->getNumEntriesInLocalRow(i * blkSize + j);
      }
    }

    // MultiJagged: Grab the target rows per process from the Heuristic to use unless the Zoltan2 list says otherwise
    if (algo == "multijagged" && !Zoltan2Params.isParameter("mj_premigration_coordinate_count")) {
      LO heuristicTargetRowsPerProcess = Get<LO>(level, "repartition: heuristic target rows per process");
      Zoltan2Params.set("mj_premigration_coordinate_count", heuristicTargetRowsPerProcess);
    }
    const bool writeZoltan2DebuggingFiles = Zoltan2Params.get("mj_debug", false);
    Zoltan2Params.remove("mj_debug");

    std::vector<int> strides;
    std::vector<const real_type*> weights(1, weightsPerRow.getRawPtr());

    RCP<const Teuchos::MpiComm<int> > dupMpiComm            = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(rowMap->getComm()->duplicate());
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > zoltanComm = dupMpiComm->getRawMpiComm();

    InputAdapterType adapter(coords, weights, strides);
    RCP<ProblemType> problem(new ProblemType(&adapter, &Zoltan2Params, (*zoltanComm)()));

    {
      SubFactoryMonitor m1(*this, "Zoltan2 " + toString(algo), level);
      if (writeZoltan2DebuggingFiles)
        adapter.generateFiles(("mj_debug.lvl_" + std::to_string(level.GetLevelID())).c_str(), *(rowMap->getComm()));
      problem->solve();
    }

    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, false);
    ArrayRCP<GO> decompEntries                         = decomposition->getDataNonConst(0);

    const typename InputAdapterType::part_t* parts = problem->getSolution().getPartListView();

    for (GO i = 0; i < numElements; i++) {
      int partNum = parts[i];

      for (LO j = 0; j < blkSize; j++)
        decompEntries[i * blkSize + j] = partNum;
    }

    Set(level, "Partition", decomposition);

  } else {
    GO numElements = rowMap->getLocalNumElements();

    typedef Zoltan2::XpetraCrsGraphAdapter<CrsGraph> InputAdapterType;
    typedef Zoltan2::PartitioningProblem<InputAdapterType> ProblemType;

    RCP<const Teuchos::MpiComm<int> > dupMpiComm            = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(rowMap->getComm()->duplicate());
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > zoltanComm = dupMpiComm->getRawMpiComm();

    InputAdapterType adapter(A->getCrsGraph());
    RCP<ProblemType> problem(new ProblemType(&adapter, &Zoltan2Params, (*zoltanComm)()));

    {
      SubFactoryMonitor m1(*this, "Zoltan2 " + toString(algo), level);
      problem->solve();
    }

    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, false);
    ArrayRCP<GO> decompEntries                         = decomposition->getDataNonConst(0);

    const typename InputAdapterType::part_t* parts = problem->getSolution().getPartListView();

    // For blkSize > 1, ignore solution for every row but the first ones in a block.
    for (GO i = 0; i < numElements / blkSize; i++) {
      int partNum = parts[i * blkSize];

      for (LO j = 0; j < blkSize; j++)
        decompEntries[i * blkSize + j] = partNum;
    }

    Set(level, "Partition", decomposition);
  }
}

}  // namespace MueLu

#endif  // if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)

#endif  // MUELU_ZOLTAN2INTERFACE_DEF_HPP
