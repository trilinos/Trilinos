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
 * MueLu_EmergencyAggregationAlgorithm_decl.hpp
 *
 *  Created on: Sep 18, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_EMERGENCYAGGREGATIONALGORITHM_DECL_HPP_
#define MUELU_EMERGENCYAGGREGATIONALGORITHM_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_AggregationAlgorithmBase.hpp"
#include "MueLu_EmergencyAggregationAlgorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"

namespace MueLu {
  /*!
    @class EmergencyAggregationAlgorithm class.
    @brief Algorithm for coarsening a graph with uncoupled aggregation.
    aggregate all leftover nodes into small "emergency" aggregates
  */

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class EmergencyAggregationAlgorithm : public MueLu::AggregationAlgorithmBase<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {
#undef MUELU_EMERGENCYAGGREGATIONALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    EmergencyAggregationAlgorithm(RCP<const FactoryBase> const &graphFact = Teuchos::null);

    //! Destructor.
    virtual ~EmergencyAggregationAlgorithm() { }

    //@}


    //! @name Aggregation methods.
    //@{

    /*! @brief Local aggregation. */

    LocalOrdinal BuildAggregates(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat, Teuchos::ArrayRCP<unsigned int> & coarse_aggStat) const;
    //@}

  }; //class MaxLinkAggregationAlgorithm

} //namespace MueLu

#define MUELU_EMERGENCYAGGREGATIONALGORITHM_SHORT
#endif /* MUELU_EMERGENCYAGGREGATIONALGORITHM_DECL_HPP_ */
