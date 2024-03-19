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
#ifndef MUELU_LWGRAPH_KOKKOS_DECL_HPP
#define MUELU_LWGRAPH_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include <Xpetra_ConfigDefs.hpp>  // global_size_t

#include "MueLu_LWGraphBase.hpp"
#include "MueLu_LWGraph_kokkos_fwd.hpp"

namespace MueLu {

/*!
  @class LWGraph_kokkos
  @brief Lightweight MueLu representation of a compressed row storage graph

  This class is lightweight in the sense that it holds to local graph
  information. These were built without using fillComplete.
 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
class LWGraph_kokkos : public MueLu::LWGraphBase<LocalOrdinal, GlobalOrdinal, Node, false> {
  using LWGraphBase<LocalOrdinal, GlobalOrdinal, Node, false>::LWGraphBase;
};

}  // namespace MueLu

#define MUELU_LWGRAPH_KOKKOS_SHORT
#endif  // MUELU_LWGRAPH_KOKKOS_DECL_HPP
