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
#ifndef MUELU_AGGREGATES_KOKKOS_DEF_HPP
#define MUELU_AGGREGATES_KOKKOS_DEF_HPP

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Utilities_kokkos_decl.hpp"
#include "MueLu_Aggregates_kokkos_decl.hpp"

namespace MueLu {

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::Aggregates_kokkos(LWGraph_kokkos graph) {
    nAggregates_  = 0;

    vertex2AggId_ = LOVectorFactory::Build(graph.GetImportMap());
    vertex2AggId_->putScalar(MUELU_UNAGGREGATED);

    procWinner_ = LOVectorFactory::Build(graph.GetImportMap());
    procWinner_->putScalar(MUELU_UNASSIGNED);

    isRoot_ = Teuchos::ArrayRCP<bool>(graph.GetImportMap()->getNodeNumElements(), false);

    // slow but safe, force TentativePFactory to build column map for P itself
    aggregatesIncludeGhosts_ = true;
  }

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::Aggregates_kokkos(const RCP<const Map>& map) {
    nAggregates_ = 0;

    vertex2AggId_ = LOVectorFactory::Build(map);
    vertex2AggId_->putScalar(MUELU_UNAGGREGATED);

    procWinner_ = LOVectorFactory::Build(map);
    procWinner_->putScalar(MUELU_UNASSIGNED);

    isRoot_ = Teuchos::ArrayRCP<bool>(map->getNodeNumElements(), false);

    // slow but safe, force TentativePFactory to build column map for P itself
    aggregatesIncludeGhosts_ = true;
  }

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::aggregates_sizes_type::const_type
  Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::ComputeAggregateSizes(bool forceRecompute, bool cacheSizes) const {
    if (aggregateSizes_.size() && !forceRecompute) {
      return aggregateSizes_;

    } else {

      // invalidate previous sizes
      aggregateSizes_ = aggregates_sizes_type("aggregates", 0);

      aggregates_sizes_type aggregateSizes("aggregates", nAggregates_);

      int myPID = vertex2AggId_->getMap()->getComm()->getRank();

      auto vertex2AggId = vertex2AggId_->template getLocalView<DeviceType>();
      auto procWinner   = procWinner_  ->template getLocalView<DeviceType>();

      typename AppendTrait<decltype(aggregateSizes_), Kokkos::Atomic>::type aggregateSizesAtomic = aggregateSizes;
      Kokkos::parallel_for("Aggregates:ComputeAggregateSizes:for", procWinner.size(), KOKKOS_LAMBDA(const LO k) {
        if (procWinner(k, 0) == myPID)
          aggregateSizesAtomic(vertex2AggId(k, 0))++;
      });

      if (cacheSizes)
        aggregateSizes_ = aggregateSizes;

      return aggregateSizes;
    }

  } //ComputeAggSizesNodes

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::local_graph_type
  Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::GetGraph() const {
    throw "Not implemented";
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  std::string Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::description() const {
    return BaseClass::description() + "{nGlobalAggregates = " + toString(GetNumGlobalAggregates()) + "}";
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::print(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Statistics0)
      out0 << "Global number of aggregates: " << GetNumGlobalAggregates() << std::endl;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::GetNumGlobalAggregates() const {
    LO nAggregates = GetNumAggregates();
    GO nGlobalAggregates;
    MueLu_sumAll(vertex2AggId_->getMap()->getComm(), (GO)nAggregates, nGlobalAggregates);
    return nGlobalAggregates;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>> >
  Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::GetMap() const {
    return vertex2AggId_->getMap();
  }

} //namespace MueLu

#endif // MUELU_AGGREGATES_KOKKOS_DEF_HPP
