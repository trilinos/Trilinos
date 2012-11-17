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
/*
 * MueLu_VariableTransferFactory_def.hpp
 *
 *  Created on: Jul 30, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AGGSTATTRANSFERFACTORY_DEF_HPP_
#define MUELU_AGGSTATTRANSFERFACTORY_DEF_HPP_

#include "MueLu_AggStatTransferFactory_decl.hpp"

#include "MueLu_CheapAggregationAlgorithm.hpp"  // needed for NodeState enum

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  AggStatTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AggStatTransferFactory(std::string const & varName)
    : varName_(varName)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  AggStatTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~AggStatTransferFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggStatTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, varName_);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggStatTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level &coarseLevel) const {

    FactoryMonitor m(*this, "AggStatTransferFactory", coarseLevel);

    // TODO find something smarter than distinction by variable name
    // for example one could use the underlaying type of the variable
    // Therefor we have to add the functionality to the Level class.
    // not sure we wanna do this. -> do decided later
    if (varName_ == "coarseAggStat") {
      Teuchos::ArrayRCP<unsigned int> data = Get<Teuchos::ArrayRCP<unsigned int> >(fineLevel,varName_);
      Set<Teuchos::ArrayRCP<unsigned int> >(coarseLevel, varName_, data);
    }

  } //Build

} // namespace MueLu


#endif /* MUELU_AGGSTATTRANSFERFACTORY_DEF_HPP_ */
