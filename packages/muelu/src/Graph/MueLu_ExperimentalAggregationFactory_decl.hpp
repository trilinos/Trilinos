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
 * MueLu_ExperimentalAggregationFactory_decl.hpp
 *
 *  Created on: Jul 26, 2012
 *      Author: wiesner
 */

#ifndef MUELU_EXPERIMENTALAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_EXPERIMENTALAGGREGATIONFACTORY_DECL_HPP_

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ExperimentalAggregationFactory_fwd.hpp"

#include "MueLu_CheapAggregationAlgorithm.hpp"
#include "MueLu_CheapAggregationAlgorithm_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class ExperimentalAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_EXPERIMENTALAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  ExperimentalAggregationFactory();

  //! Destructor.
  virtual ~ExperimentalAggregationFactory() { }

  //@}

  //! @name Set/get methods.
  //@{

  // Options algo1
  void SetOrdering(Ordering ordering) { algo1_->SetOrdering(ordering); }
  void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) { algo1_->SetMaxNeighAlreadySelected(maxNeighAlreadySelected); }
  AggOptions::Ordering GetOrdering() const { return algo1_->GetOrdering(); }
  int GetMaxNeighAlreadySelected() const { return algo1_->GetMaxNeighAlreadySelected(); }

  // Options shared algo1 and algo2
  void SetMinNodesPerAggregate(int minNodesPerAggregate) { algo1_->SetMinNodesPerAggregate(minNodesPerAggregate); }
  int GetMinNodesPerAggregate() const { return algo1_->GetMinNodesPerAggregate(); }

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

  //! Algorithms
  RCP<MueLu::CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > algo1_;

}; // class ExperimentalAggregationFactory

}

#define MUELU_EXPERIMENTALAGGREGATIONFACTORY_SHORT
#endif /* MUELU_EXPERIMENTALAGGREGATIONFACTORY_DECL_HPP_ */
