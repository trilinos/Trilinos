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
#ifndef MUELU_AGGREGATES_DEF_HPP
#define MUELU_AGGREGATES_DEF_HPP

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_BlockedMultiVector.hpp>
#include <Xpetra_BlockedVector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Aggregates_decl.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Utilities_decl.hpp" // MueLu_sumAll

namespace MueLu {

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Aggregates<LocalOrdinal, GlobalOrdinal, Node>::Aggregates(const GraphBase & graph) {
    nAggregates_  = 0;

    vertex2AggId_ = LOMultiVectorFactory::Build(graph.GetImportMap(), 1);
    vertex2AggId_->putScalar(MUELU_UNAGGREGATED);

    procWinner_ = LOVectorFactory::Build(graph.GetImportMap());
    procWinner_->putScalar(MUELU_UNASSIGNED);

    isRoot_ = Teuchos::ArrayRCP<bool>(graph.GetImportMap()->getNodeNumElements(), false);

    // slow but safe, force TentativePFactory to build column map for P itself
    aggregatesIncludeGhosts_ = true;
  }

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Aggregates<LocalOrdinal, GlobalOrdinal, Node>::Aggregates(const RCP<const Map> & map) {
    nAggregates_ = 0;

    vertex2AggId_ = LOMultiVectorFactory::Build(map, 1);
    vertex2AggId_->putScalar(MUELU_UNAGGREGATED);

    procWinner_ = LOVectorFactory::Build(map);
    procWinner_->putScalar(MUELU_UNASSIGNED);

    isRoot_ = Teuchos::ArrayRCP<bool>(map->getNodeNumElements(), false);

    // slow but safe, force TentativePFactory to build column map for P itself
    aggregatesIncludeGhosts_ = true;
  }

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<LocalOrdinal>  Aggregates<LocalOrdinal, GlobalOrdinal, Node>::ComputeAggregateSizes(bool forceRecompute) const {

    if (aggregateSizes_ != Teuchos::null && !forceRecompute) {

      return aggregateSizes_;

    } else {

      //invalidate previous sizes.
      aggregateSizes_ = Teuchos::null;

      Teuchos::ArrayRCP<LO> aggregateSizes;
      aggregateSizes = Teuchos::ArrayRCP<LO>(nAggregates_,0);
      int myPid = vertex2AggId_->getMap()->getComm()->getRank();
      Teuchos::ArrayRCP<LO> procWinner   = procWinner_->getDataNonConst(0);
      Teuchos::ArrayRCP<LO> vertex2AggId = vertex2AggId_->getDataNonConst(0);
      LO size = procWinner.size();

      //for (LO i = 0; i < nAggregates_; ++i) aggregateSizes[i] = 0;
      for (LO k = 0; k < size; ++k ) {
        if (procWinner[k] == myPid) aggregateSizes[vertex2AggId[k]]++;
      }

      aggregateSizes_ = aggregateSizes;

      return aggregateSizes;
    }

  } //ComputeAggSizesNodes

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string Aggregates<LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;
    out << BaseClass::description();
    out << "{nGlobalAggregates = " << GetNumGlobalAggregates() << "}";
    return out.str();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Aggregates<LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Statistics1)
      out0 << "Global number of aggregates: " << GetNumGlobalAggregates() << std::endl;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal Aggregates<LocalOrdinal, GlobalOrdinal, Node>::GetNumGlobalAggregates() const {
    LO nAggregates = GetNumAggregates();
    GO nGlobalAggregates; MueLu_sumAll(vertex2AggId_->getMap()->getComm(), (GO)nAggregates, nGlobalAggregates);
    return nGlobalAggregates;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> > Aggregates<LocalOrdinal, GlobalOrdinal, Node>::GetMap() const {
    return vertex2AggId_->getMap();
  }

} //namespace MueLu

#endif // MUELU_AGGREGATES_DEF_HPP
