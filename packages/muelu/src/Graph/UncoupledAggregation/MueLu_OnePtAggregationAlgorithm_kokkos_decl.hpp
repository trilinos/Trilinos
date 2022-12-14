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
#ifndef MUELU_ONEPTAGGREGATIONALGORITHM_KOKKOS_DECL_HPP
#define MUELU_ONEPTAGGREGATIONALGORITHM_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include "MueLu_OnePtAggregationAlgorithm_kokkos_fwd.hpp"

#include "MueLu_Aggregates_kokkos_fwd.hpp"
#include "MueLu_AggregationAlgorithmBase_kokkos.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_LWGraph_kokkos.hpp"

namespace MueLu {
  /*!
    @class OnePtAggregationAlgorithm class.
    @brief Algorithm for coarsening a graph with uncoupled aggregation.
    keep special marked nodes as singleton node aggregates over all multigrid levels

    @ingroup Aggregation

    ### Idea ###
    The user can mark some nodes as ONEPT to build some single node aggregates.
    This can be very useful for certain applications. We build single node aggregates
    for nodes with the state ONEPT. Then, the state is changed to ignored.
    The OnePtAggregationAlgorithm should run before the Phase1AggregationAlgorithm.

    ### Comments ###
    Only nodes with state ONEPT are changed to IGNORED.

  */

  template<class LocalOrdinal = DefaultLocalOrdinal,
           class GlobalOrdinal = DefaultGlobalOrdinal,
           class Node = DefaultNode>
  class OnePtAggregationAlgorithm_kokkos :
    public MueLu::AggregationAlgorithmBase_kokkos<LocalOrdinal,GlobalOrdinal,Node> {
#undef MUELU_ONEPTAGGREGATIONALGORITHM_KOKKOS_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:
    using device_type  = typename LWGraph_kokkos::device_type;
    using memory_space = typename LWGraph_kokkos::memory_space;
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    OnePtAggregationAlgorithm_kokkos(RCP<const FactoryBase> const &graphFact = Teuchos::null);

    //! Destructor.
    virtual ~OnePtAggregationAlgorithm_kokkos() { }

    //@}


    //! @name Aggregation methods.
    //@{

    /*! @brief Local aggregation. */

    void BuildAggregates(Teuchos::ParameterList const & params,
                         LWGraph_kokkos const & graph,
                         Aggregates_kokkos & aggregates,
                         Kokkos::View<unsigned*, device_type>& aggStat,
                         LO& numNonAggregatedNodes) const;
    //@}


  }; //class OnePtAggregationAlgorithm_kokkos

} //namespace MueLu

#define MUELU_ONEPTAGGREGATIONALGORITHM_KOKKOS_SHORT
#endif // MUELU_ONEPTAGGREGATIONALGORITHM_KOKKOS_DECL_HPP
