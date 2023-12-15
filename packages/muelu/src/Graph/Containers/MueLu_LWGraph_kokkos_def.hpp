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
#ifndef MUELU_LWGRAPH_KOKKOS_DEF_HPP
#define MUELU_LWGRAPH_KOKKOS_DEF_HPP

#include <Kokkos_Core.hpp>

#include <Teuchos_ArrayView.hpp>
#include <Xpetra_Map.hpp>

#include "MueLu_LWGraph_kokkos_decl.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
void LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>>::
    print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
  if (verbLevel & Debug) {
    auto graph             = lclLWGraph_.getGraph();
    RCP<const Map> col_map = importMap_.is_null() ? domainMap_ : importMap_;
    int mypid              = col_map->getComm()->getRank();

    {
      std::ostringstream ss;
      ss << "[pid " << mypid << "] num entries=" << graph.entries.size();
      out << ss.str() << std::endl;
    }

    const size_t numRows = graph.numRows();
    auto rowPtrs         = graph.row_map;
    auto columns         = graph.entries;
    for (size_t i = 0; i < numRows; ++i) {
      std::ostringstream ss;
      ss << "[pid " << mypid << "] row " << domainMap_->getGlobalElement(i) << ":";
      ss << " (numEntries=" << rowPtrs(i + 1) - rowPtrs(i) << ")";

      auto rowView = graph.rowConst(i);
      for (LO j = 0; j < rowView.length; j++) {
        ss << " " << col_map->getGlobalElement(rowView.colidx(j));
      }
      out << ss.str() << std::endl;
    }
  }
}

}  // namespace MueLu

#endif  // MUELU_LWGRAPH_KOKKOS_DEF_HPP
