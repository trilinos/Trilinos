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
#ifndef PACKAGES_MUELU_SRC_REBALANCING_MUELU_REPARTITIONHEURISTICFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_REBALANCING_MUELU_REPARTITIONHEURISTICFACTORY_DEF_HPP_


#include <algorithm>
#include <iostream>
#include <sstream>

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>

//#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_Utilities.hpp"

#include "MueLu_RAPFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_RepartitionHeuristicFactory.hpp"

namespace MueLu {

 template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
 RCP<const ParameterList> RepartitionHeuristicFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("repartition: start level");
    SET_VALID_ENTRY("repartition: min rows per proc");
    SET_VALID_ENTRY("repartition: target rows per proc");
    SET_VALID_ENTRY("repartition: max imbalance");
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",         Teuchos::null, "Factory of the matrix A");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RepartitionHeuristicFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RepartitionHeuristicFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    const Teuchos::ParameterList & pL = GetParameterList();
    // Access parameters here to make sure that we set the parameter entry flag to "used" even in case of short-circuit evaluation.
    // TODO (JG): I don't really know if we want to do this.
    const int    startLevel           = pL.get<int>   ("repartition: start level");
    const LO     minRowsPerProcess    = pL.get<LO>    ("repartition: min rows per proc");
          LO     targetRowsPerProcess = pL.get<LO>    ("repartition: target rows per proc");
    const double nonzeroImbalance     = pL.get<double>("repartition: max imbalance");

    if (targetRowsPerProcess == 0)
      targetRowsPerProcess = minRowsPerProcess;

    RCP<const FactoryBase> Afact = GetFactory("A");
    if(!Afact.is_null() && Teuchos::rcp_dynamic_cast<const RAPFactory>(Afact) == Teuchos::null &&
       Teuchos::rcp_dynamic_cast<const BlockedRAPFactory>(Afact) == Teuchos::null &&
       Teuchos::rcp_dynamic_cast<const SubBlockAFactory>(Afact) == Teuchos::null) {
      GetOStream(Warnings) <<
        "MueLu::RepartitionHeuristicFactory::Build: The generation factory for A must " \
        "be a RAPFactory or a SubBlockAFactory providing the non-rebalanced matrix information! " \
        "It specifically must not be of type Rebalance(Blocked)AcFactory or similar. " \
        "Please check the input. Make also sure that \"number of partitions\" is provided to " \
        "the Interface class and the RepartitionFactory instance.  Instead, we have a "<<Afact->description() << std::endl;
    }
    // TODO: We only need a CrsGraph. This class does not have to be templated on Scalar types.
    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");

    // ======================================================================================================
    // Determine whether partitioning is needed
    // ======================================================================================================
    // NOTE: most tests include some global communication, which is why we currently only do tests until we make
    // a decision on whether to repartition. However, there is value in knowing how "close" we are to having to
    // rebalance an operator. So, it would probably be beneficial to do and report *all* tests.

    // Test1: skip repartitioning if current level is less than the specified minimum level for repartitioning
    if (currentLevel.GetLevelID() < startLevel) {
      GetOStream(Statistics1) << "Repartitioning?  NO:" <<
          "\n  current level = " << Teuchos::toString(currentLevel.GetLevelID()) <<
          ", first level where repartitioning can happen is " + Teuchos::toString(startLevel) << std::endl;

      // a negative number of processors means: no repartitioning
      Set(currentLevel, "number of partitions", -1);

      return;
    }

    RCP<const Map> rowMap = A->getRowMap();

    RCP<const Teuchos::Comm<int> > origComm = rowMap->getComm();
    RCP<const Teuchos::Comm<int> > comm     = origComm;

    // Test 2: check whether A is actually distributed, i.e. more than one processor owns part of A
    // TODO: this global communication can be avoided if we store the information with the matrix (it is known when matrix is created)
    // TODO: further improvements could be achieved when we use subcommunicator for the active set. Then we only need to check its size
    {
      int numActiveProcesses = 0;
      MueLu_sumAll(comm, Teuchos::as<int>((A->getNodeNumRows() > 0) ? 1 : 0), numActiveProcesses);

      if (numActiveProcesses == 1) {
        GetOStream(Statistics1) << "Repartitioning?  NO:" <<
            "\n  # processes with rows = " << Teuchos::toString(numActiveProcesses) << std::endl;

        Set(currentLevel, "number of partitions", 1);
        return;
      }
    }

    bool test3 = false, test4 = false;
    std::string msg3, msg4;

    // Test3: check whether number of rows on any processor satisfies the minimum number of rows requirement
    // NOTE: Test2 ensures that repartitionning is not done when there is only one processor (it may or may not satisfy Test3)
    if (minRowsPerProcess > 0) {
      LO numMyRows = Teuchos::as<LO>(A->getNodeNumRows()), minNumRows, LOMAX = Teuchos::OrdinalTraits<LO>::max();
      LO haveFewRows = (numMyRows < minRowsPerProcess ? 1 : 0), numWithFewRows = 0;
      MueLu_sumAll(comm, haveFewRows, numWithFewRows);
      MueLu_minAll(comm, (numMyRows > 0 ? numMyRows : LOMAX), minNumRows);

      // TODO: we could change it to repartition only if the number of processors with numRows < minNumRows is larger than some
      // percentage of the total number. This way, we won't repartition if 2 out of 1000 processors don't have enough elements.
      // I'm thinking maybe 20% threshold. To implement, simply add " && numWithFewRows < .2*numProcs" to the if statement.
      if (numWithFewRows > 0)
        test3 = true;

      msg3 = "\n  min # rows per proc = " + Teuchos::toString(minNumRows) + ", min allowable = " + Teuchos::toString(minRowsPerProcess);
    }

    // Test4: check whether the balance in the number of nonzeros per processor is greater than threshold
    if (!test3) {
      GO minNnz, maxNnz, numMyNnz = Teuchos::as<GO>(A->getNodeNumEntries());
      MueLu_maxAll(comm, numMyNnz,                           maxNnz);
      MueLu_minAll(comm, (numMyNnz > 0 ? numMyNnz : maxNnz), minNnz); // min nnz over all active processors
      double imbalance = Teuchos::as<double>(maxNnz)/minNnz;

      if (imbalance > nonzeroImbalance)
        test4 = true;

      msg4 = "\n  nonzero imbalance = " + Teuchos::toString(imbalance) + ", max allowable = " + Teuchos::toString(nonzeroImbalance);
    }

    if (!test3 && !test4) {
      GetOStream(Statistics1) << "Repartitioning?  NO:" << msg3 + msg4 << std::endl;

      // A negative number of partitions means: no repartitioning
      Set(currentLevel, "number of partitions", -1);
      return;
    }

    GetOStream(Statistics1) << "Repartitioning? YES:" << msg3 + msg4 << std::endl;

    // ======================================================================================================
    // Calculate number of partitions
    // ======================================================================================================
    // FIXME Quick way to figure out how many partitions there should be (same algorithm as ML)
    // FIXME Should take into account nnz? Perhaps only when user is using min #nnz per row threshold.

    // The number of partitions is calculated by the RepartitionFactory and stored in "number of partitions" variable on
    // the current level. If this variable is already set (e.g., by another instance of RepartitionFactory) then this number
    // is used. The "number of partitions" variable serves as basic communication between the RepartitionFactory (which
    // requests a certain number of partitions) and the *Interface classes which call the underlying partitioning algorithms
    // and produce the "Partition" array with the requested number of partitions.
    int numPartitions;
    if (Teuchos::as<GO>(A->getGlobalNumRows()) < targetRowsPerProcess) {
      // System is too small, migrate it to a single processor
      numPartitions = 1;

    } else {
      // Make sure that each processor has approximately targetRowsPerProcess
      numPartitions = A->getGlobalNumRows() / targetRowsPerProcess;
    }
    numPartitions = std::min(numPartitions, comm->getSize());

    Set(currentLevel, "number of partitions", numPartitions);

    GetOStream(Statistics1) << "Number of partitions to use = " << numPartitions << std::endl;
  } // Build
} // namespace MueLu

#endif //ifdef HAVE_MPI
#endif /* PACKAGES_MUELU_SRC_REBALANCING_MUELU_REPARTITIONHEURISTICFACTORY_DEF_HPP_ */
