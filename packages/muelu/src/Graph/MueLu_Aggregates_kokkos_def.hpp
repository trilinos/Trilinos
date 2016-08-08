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
  template<class ProcWinnerType, class Vertex2AggIdType, class AggregateSizesType, class LO>
  class ComputeAggregateSizesFunctor {
  private:
    ProcWinnerType      procWinner;
    Vertex2AggIdType    vertex2AggId;
    int                 myPID;
    AggregateSizesType  aggregateSizes;

  public:
    ComputeAggregateSizesFunctor(ProcWinnerType procWinner_, Vertex2AggIdType vertex2AggId_, int myPID_, AggregateSizesType aggregateSizes_) :
      procWinner(procWinner_),
      vertex2AggId(vertex2AggId_),
      myPID(myPID_),
      aggregateSizes(aggregateSizes_)
    { }

    KOKKOS_INLINE_FUNCTION
    void operator()(const LO k) const {
      if (procWinner(k, 0) == myPID)
        aggregateSizes(vertex2AggId(k, 0))++;
    }
  };
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::aggregates_sizes_type::const_type
  Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::ComputeAggregateSizes(bool forceRecompute, bool cacheSizes) const {
    if (aggregateSizes_.size() && !forceRecompute) {
      return aggregateSizes_;

    } else {

      // invalidate previous sizes
      aggregateSizes_ = aggregates_sizes_type("aggregates", 0);

      aggregates_sizes_type aggregateSizes("aggregates", nAggregates_);

      int myPID = GetMap()->getComm()->getRank();

      auto vertex2AggId = vertex2AggId_->template getLocalView<DeviceType>();
      auto procWinner   = procWinner_  ->template getLocalView<DeviceType>();

      typename AppendTrait<decltype(aggregateSizes_), Kokkos::Atomic>::type aggregateSizesAtomic = aggregateSizes;

      ComputeAggregateSizesFunctor<decltype(procWinner), decltype(vertex2AggId), decltype(aggregateSizesAtomic), LO>
        computeAggSizesFunctor(procWinner, vertex2AggId, myPID, aggregateSizesAtomic);
      Kokkos::parallel_for("MueLu:Aggregates:ComputeAggregateSizes:for", procWinner.size(), computeAggSizesFunctor);

      if (cacheSizes)
        aggregateSizes_ = aggregateSizes;

      return aggregateSizes;
    }

  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::local_graph_type
  Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::GetGraph() const {
    typedef typename local_graph_type::row_map_type row_map_type;
    typedef typename local_graph_type::entries_type entries_type;

    int myPID = GetMap()->getComm()->getRank();

    ArrayRCP<LO> vertex2AggId = vertex2AggId_->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = procWinner_->getDataNonConst(0);

    typename aggregates_sizes_type::const_type sizes = ComputeAggregateSizes();

    int numAggregates = nAggregates_;

    typename row_map_type::non_const_type rows("row_map", numAggregates+1);       // rows(0) = 0 automatically
    for (LO i = 0; i < nAggregates_; i++) // TODO: replace by parallel_scan
      rows(i+1) = rows(i) + sizes(i);

    aggregates_sizes_type offsets("offsets", numAggregates);
    for (LO i = 0; i < numAggregates; i++) // TODO: replace by parallel_for
      offsets(i) = rows(i);

    typename entries_type::non_const_type cols("entries", rows(nAggregates_));
    for (LO i = 0; i < procWinner.size(); i++)
      if (procWinner[i] == myPID)
        cols(offsets(vertex2AggId[i])++) = i;

    return local_graph_type(cols, rows);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  std::string Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::description() const {
    return BaseClass::description() + "{nGlobalAggregates = " + toString(GetNumGlobalAggregates()) + "}";
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::print(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Statistics1)
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
