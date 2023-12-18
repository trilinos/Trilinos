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
#ifndef MUELU_LOCALLWGRAPH_KOKKOS_DEF_HPP
#define MUELU_LOCALLWGRAPH_KOKKOS_DEF_HPP

#include <Kokkos_Core.hpp>

#include <Teuchos_ArrayView.hpp>
#include <Xpetra_Map.hpp>

#include "MueLu_LocalLWGraph_kokkos_decl.hpp"

namespace MueLu {

namespace {  // anonymous

template <class LocalOrdinal, class RowType>
class MaxNumRowEntriesFunctor {
 public:
  MaxNumRowEntriesFunctor(RowType rowPointers)
    : rowPointers_(rowPointers) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const LocalOrdinal i, size_t& maxLength) const {
    size_t d = rowPointers_(i + 1) - rowPointers_(i);

    maxLength = (d > maxLength ? d : maxLength);
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile size_t& dest, const volatile size_t& src) {
    dest = (dest > src ? dest : src);
  }

  KOKKOS_INLINE_FUNCTION
  void init(size_t& initValue) {
    initValue = 0;
  }

 private:
  RowType rowPointers_;
};

}  // namespace

template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
LocalLWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>>::
    LocalLWGraph_kokkos(const local_graph_type& graph,
                        const RCP<const map_type>& domainMap)
  : graph_(graph) {
  minLocalIndex_ = domainMap->getMinLocalIndex();
  maxLocalIndex_ = domainMap->getMaxLocalIndex();

  MaxNumRowEntriesFunctor<LO, typename local_graph_type::row_map_type> maxNumRowEntriesFunctor(graph_.row_map);
  Kokkos::parallel_reduce("MueLu:LocalLWGraph:LWGraph:maxnonzeros", range_type(0, graph_.numRows()), maxNumRowEntriesFunctor, maxNumRowEntries_);
}

}  // namespace MueLu

#endif  // MUELU_LWGRAPH_KOKKOS_DEF_HPP
