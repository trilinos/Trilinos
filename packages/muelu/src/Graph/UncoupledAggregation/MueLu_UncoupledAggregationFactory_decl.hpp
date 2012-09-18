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

 * MueLu_UncoupledAggregationFactory_decl.hpp
 *
 *  Created on: Sep 17, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_UNCOUPLEDAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_DECL_HPP_


#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"

//#include "MueLu_CheapAggregationAlgorithm.hpp"
//#include "MueLu_CheapAggregationAlgorithm_fwd.hpp"


#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_OnePtAggregationAlgorithm_fwd.hpp"
#include "MueLu_UncoupledAggregationAlgorithm_fwd.hpp"
#include "MueLu_MaxLinkAggregationAlgorithm_fwd.hpp"
#include "MueLu_EmergencyAggregationAlgorithm_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class UncoupledAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_UNCOUPLEDAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  UncoupledAggregationFactory(RCP<const FactoryBase> graphFact = Teuchos::null);

  //! Destructor.
  virtual ~UncoupledAggregationFactory() { }

  //@}

  //! @name Set/get methods.
  //@{

  // Options algo1
  void SetOrdering(AggOptions::UncoupledOrdering ordering) {
    for(size_t a = 0; a < algos_.size(); a++) {
      algos_[a]->SetOrdering(ordering);
    }
  }
  void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) {
    for(size_t a = 0; a < algos_.size(); a++) {
      algos_[a]->SetMaxNeighAlreadySelected(maxNeighAlreadySelected);
    }
  }
  AggOptions::UncoupledOrdering GetOrdering() const {
    return algos_[0]->GetOrdering();
  }
  int GetMaxNeighAlreadySelected() const {
    return algos_[0]->GetMaxNeighAlreadySelected();
  }

  // Options shared algo1 and algo2
  void SetMinNodesPerAggregate(int minNodesPerAggregate) {
    for(size_t a = 0; a < algos_.size(); a++) {
      algos_[a]->SetMinNodesPerAggregate(minNodesPerAggregate);
    }
  }
  int GetMinNodesPerAggregate() const {
    return algos_[0]->GetMinNodesPerAggregate();
  }

  //@}

  //! Input
  //@{

  void DeclareInput(Level &currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates. */
  void Build(Level &currentLevel) const;

  //@}

private:

  //! Graph Factory
  RCP<const FactoryBase> graphFact_;

  //! Algorithms
  std::vector<RCP<MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > > algos_;

}; // class UncoupledAggregationFactory

}

#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_SHORT
#endif /* MUELU_UNCOUPLEDAGGREGATIONFACTORY_DECL_HPP_ */
