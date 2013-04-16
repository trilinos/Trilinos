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
#include <set>

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

    validParamList->set< RCP<const FactoryBase> >("A",                      Teuchos::null, "Factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",            Teuchos::null, "Factory of the coordinates");
    validParamList->set< RCP<const FactoryBase> >("number of partitions",   Teuchos::null, "(advanced) Factory computing the number of partitions");
    validParamList->set< int >                   ("rowWeight",                          0, "Default weight to rows (total weight = nnz + rowWeight");
    validParamList->set< std::string >           ("algorithm",              "multijagged", "Zoltan2 partitioning algorithm (multijagged,rcb)");
    validParamList->set< std::string >           ("pqparts",                       "2k3m", "Algorithm for deciding the pqparts split");

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

    TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->lib() == Xpetra::UseEpetra, Exceptions::RuntimeError, "Zoltan2 does not work with Epetra at the moment. Please use Zoltan through ZoltanInterface");

    if (!IsAvailable(level, "Coordinates")) {
      std::cout << GetFactory("Coordinates") << std::endl;

      level.print(*getOStream());
      throw Exceptions::RuntimeError("MueLu::ZoltanInterface::Build(): no coordinates available");
    }

    RCP<MultiVector> coords = Get< RCP<MultiVector> >(level, "Coordinates");
    size_t              dim = coords->getNumVectors();

    // Check that the number of local coordinates is consistent with the #rows in A
    std::string msg = "MueLu::Zoltan2Interface::Build : coordinate vector length is incompatible with number of rows in A.  The vector length should be the same as the number of mesh points.";
    TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->getNodeNumElements()/blkSize != coords->getLocalLength(), Exceptions::Incompatible, msg);

    RCP<const Map> map = coords->getMap();

    if (numParts == 1) {
      // Single processor, decomposition is trivial: all zeros
      RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), true);
      Set(level, "Partition", decomposition);
      return;
    }

    GO numElements = map->getNodeNumElements();
    std::vector<const SC*> values(dim), weights(1);
    std::vector<int>       strides;

    for (size_t k = 0; k < dim; k++)
      values[k] = coords->getData(k).get();

    const ParameterList& pL = GetParameterList();
    int rowWeight = pL.get<int>("rowWeight");
    GetOStream(Runtime0,0) << "Using weights formula: nnz + " << rowWeight << std::endl;

    Array<SC> weightsPerRow(numElements);
    for (LO i = 0; i < numElements; i++)
      for (LO j = 0; j < blkSize; j++) {
        weightsPerRow[i] += A->getNumEntriesInLocalRow(i*blkSize+j);
        // Zoltan2 pqJagged gets as good partitioning as Zoltan RCB in terms of nnz
        // but Zoltan also gets a good partioning in rows, which sometimes does not
        // happen for Zoltan2. So here is an attempt to get a better row partitioning
        // without significantly screwing up nnz partitioning
        // NOTE: no good heuristic here, the value was chosen almost randomly
        weightsPerRow[i] += rowWeight;
      }
    weights[0] = weightsPerRow.getRawPtr();

    Teuchos::ParameterList params;

    std::string algo = pL.get<std::string>("algorithm");
    TEUCHOS_TEST_FOR_EXCEPTION(algo != "multijagged" && algo != "rcb", Exceptions::RuntimeError, "Unknown partitioning algorithm: \"" << algo << "\"");
    if (algo == "multijagged") {
      params.set("algorithm",             "multijagged");
      std::string pqalgo = pL.get<std::string>("pqparts");
      params.set("pqParts",               getPQParts(pqalgo, numParts, dim));

    } else if (algo == "rcb") {
      params.set("algorithm",             "rcb");
      params.set("num_global_parts",      numParts);
    }

    params.set("partitioning_approach",   "partition");
    // NOTE: "compute metrics" is buggy in Zoltan2, especially for large runs
    // params.set("compute_metrics",         "true");

    GetOStream(Runtime0,0) << "Zoltan2 parameters:" << std::endl << "----------" << std::endl << params << "----------" << std::endl;

    typedef Zoltan2::BasicCoordinateInput<Zoltan2::BasicUserTypes<SC,GO,LO,GO> > InputAdapterType;
    typedef Zoltan2::PartitioningProblem<InputAdapterType> ProblemType;

    InputAdapterType adapter(numElements, map->getNodeElementList().getRawPtr(), values, strides, weights, strides);

    const Teuchos::MpiComm<int>& comm = static_cast<const Teuchos::MpiComm<int>& >(*map->getComm());
    RCP<ProblemType> problem(new ProblemType(&adapter, &params, *comm.getRawMpiComm()));

    problem->solve();

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(A->getRowMap(), false);
    ArrayRCP<GO>                      decompEntries = decomposition->getDataNonConst(0);

    const zoltan2_partId_t * parts = problem->getSolution().getPartList();

    // KDDKDD NEW:  At present, Zoltan2 does not guarantee that the
    // KDDKDD NEW:  parts in getPartList() are listed in the same order
    // KDDKDD NEW:  as the input.  Using getIdList() compensates for
    // KDDKDD NEW:  differences in the order.  Eventually, Zoltan2 will
    // KDDKDD NEW:  guarantee identical ordering; at that time, all code
    // KDDKDD NEW:  marked with "KDDKDD NEW" can be reverted to the code
    // KDDKDD NEW:  marked with "KDDKDD OLD".
    const GO * zgids = problem->getSolution().getIdList();  // KDDKDD  NEW

    for (GO i = 0; i < numElements; i++) {
      //GO localID = i;   // KDDKDD OLD
      GO localID = map->getLocalElement(zgids[i]); // KDDKDD NEW
      int partNum = parts[i];

      for (LO j = 0; j < blkSize; j++)
        decompEntries[localID*blkSize + j] = partNum;
    }

    Set(level, "Partition", decomposition);

  } //Build()

  // Algorithm is kind of complex, but the details are not important
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string Zoltan2Interface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::getPQParts(const std::string& algo, GO numParts, size_t dim) const {
    std::ostringstream pqParts;

    TEUCHOS_TEST_FOR_EXCEPTION(algo != "2k3m" && algo != "2k", Exceptions::RuntimeError, "Unknown pqparts determination algorithm: \"" << algo << "\"");

    if (dim == 1) {
      pqParts << numParts;
      return pqParts.str();
    }

    // We would like to find a good recursive representation of the numProcs. It is not
    // possible for all values, so we try to find a close value which has one
    if (algo == "2k3m") {
      // Here is one of the possible algorithms: try to find a closest number
      // of form 2^k*3^m which is smaller than numParts. Generally, this makes the number
      // of processors off by less than 15%
      int i2 = -1, m2 = Teuchos::as<int>(floor(log(numParts)/log(2)));
      int i3 = -1, m3 = Teuchos::as<int>(floor(log(numParts)/log(3)));
      int d = Teuchos::as<int>(1e+9);

      for (int i = 0; i <= m2; i++)
        for (int j = 0; j <= m3; j++) {
          int k = Teuchos::as<int>(std::pow(2.,i) * std::pow(3.,j));
          if (k <= numParts && (numParts - k < d)) {
            d = numParts - k;
            i2 = i;
            i3 = j;
          }
        }
      if (d) {
        numParts -= d;
        GetOStream(Runtime0,0) << "Changing number of partitions to " << numParts << " in order to improve Zoltan2 pqJagged algorithm" << std::endl;
      }

      std::multiset<int> vals;
      std::multiset<int>::const_reverse_iterator rit, rit1;

      // Step 1: initialization with 1,2,3 factors
      for (int i = 1; i <= i2; i++)           vals.insert(2);
      for (int i = 1; i <= i3; i++)           vals.insert(3);
      while (vals.size() < dim)               vals.insert(1);

      // Step 2: create factor products in order to limit ourselves to 2*dim factors
      while (vals.size() > 2*dim) {
        int v1 = *vals.begin(); vals.erase(vals.begin());
        int v2 = *vals.begin(); vals.erase(vals.begin());
        vals.insert(v1*v2);
      }
      // Save the current state for the possible fallback
      std::multiset<int> vals_copy = vals;

      // Step 3: try to do some further aggregation
      const int threshold_factor = 4;
      for (size_t k = 0; dim+k < 5 && vals.size() > dim; k++) {
        int v1 = *vals.begin(), v2 = *(++vals.begin());
        rit = vals.rbegin()++;
        if (k+1 != dim)
          rit++;
        if (v1*v2 <= threshold_factor * (*rit)) {
          vals.erase(vals.begin());
          vals.erase(vals.begin());
          vals.insert(v1*v2);
        }
      }

      // Step 4: check if the aggregation is good enough
      if (vals.size() > dim && vals.size() < 2*dim) {
        rit = vals.rbegin();
        for (size_t k = 0; k < dim; k++, rit++);
        for (; rit != vals.rend() && (*rit <= 6); rit++);
        if (rit != vals.rend()) {
          // We don't like this factorization, fallback
          vals = vals_copy;
        }
      }

      // Step 5: save factors in a string
      // Tricky: first <dim> factors are printed in reverse order
      rit1 = vals.rbegin();
      for (size_t k = 0; k < dim; k++) {
        rit = vals.rbegin();
        for (size_t j = 1; j+k < dim; j++)
          rit++;
        if (k)
          pqParts << ",";
        pqParts << *rit;
        rit1++;
      }
      for (; rit1 != vals.rend(); rit1++)
        pqParts << "," << *rit1;

    } else if (algo == "2k") {
      // Here is one of the possible algorithms: try to find a closest number
      // of form 2^k which is smaller than numParts. This makes the number
      // of processors off by less than 50%
      int i2 = Teuchos::as<int>(floor(log(numParts)/log(2)));

      int d = numParts - Teuchos::as<int>(std::pow(2.,i2));
      if (d) {
        numParts -= d;
        GetOStream(Runtime0,0) << "Changing number of partitions to " << numParts << " in order to improve Zoltan2 pqJagged algorithm" << std::endl;
      }

      std::multiset<int> vals;

      // Step 1: initialization with 1,2 factors
      for (int i = 1; i <= i2; i++)           vals.insert(2);
      while (vals.size() < dim)               vals.insert(1);

      // Step 2: save factors in a string
      pqParts << *(vals.begin());
      for (std::multiset<int>::const_iterator it = (++vals.begin()); it != vals.end(); it++)
        pqParts << "," << *it;
    }

    return pqParts.str();
  }

} //namespace MueLu

#endif //if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)

#endif // MUELU_ZOLTAN2INTERFACE_DEF_HPP
