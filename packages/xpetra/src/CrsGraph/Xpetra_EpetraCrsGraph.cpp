// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
#include "Xpetra_EpetraCrsGraph.hpp"

#include "Xpetra_Exceptions.hpp"
#include "Xpetra_Utils.hpp"
#include "Xpetra_EpetraExport.hpp"
#include "Xpetra_EpetraImport.hpp"

namespace Xpetra {

  // TODO: move that elsewhere
  template<class GlobalOrdinal, class Node>
  const Epetra_CrsGraph & toEpetra(const RCP< const CrsGraph<int, GlobalOrdinal, Node> > &graph) {
    XPETRA_RCP_DYNAMIC_CAST(const EpetraCrsGraphT<GlobalOrdinal XPETRA_COMMA Node>, graph, epetraGraph, "toEpetra");
    return *(epetraGraph->getEpetra_CrsGraph());
  }

  // TODO: move that elsewhere
  template<class GlobalOrdinal, class Node>
  RCP<const CrsGraph<int, GlobalOrdinal, Node> >
  toXpetra (const Epetra_CrsGraph &g)
  {
    RCP<const Epetra_CrsGraph> const_graph = rcp (new Epetra_CrsGraph (g));
    RCP<Epetra_CrsGraph> graph =
      Teuchos::rcp_const_cast<Epetra_CrsGraph> (const_graph);
    return rcp (new Xpetra::EpetraCrsGraphT<GlobalOrdinal,Node>(graph));
  }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
    (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
  template class EpetraCrsGraphT<int, Xpetra::EpetraNode >;
  template RCP< const CrsGraph<int, int, Xpetra::EpetraNode > > toXpetra<int, Xpetra::EpetraNode>(const Epetra_CrsGraph &g);
  template const Epetra_CrsGraph & toEpetra<int, Xpetra::EpetraNode >(const RCP< const CrsGraph<int, int, Xpetra::EpetraNode > > &graph);
#endif
#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraCrsGraphT<int, Kokkos::Compat::KokkosSerialWrapperNode >;
template RCP< const CrsGraph<int, int, Kokkos::Compat::KokkosSerialWrapperNode > > toXpetra<int, Kokkos::Compat::KokkosSerialWrapperNode>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<int, Kokkos::Compat::KokkosSerialWrapperNode >(const RCP< const CrsGraph<int, int, Kokkos::Compat::KokkosSerialWrapperNode > > &graph);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraCrsGraphT<int, Kokkos::Compat::KokkosThreadsWrapperNode>;
template RCP< const CrsGraph<int, int, Kokkos::Compat::KokkosThreadsWrapperNode > > toXpetra<int, Kokkos::Compat::KokkosThreadsWrapperNode>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<int, Kokkos::Compat::KokkosThreadsWrapperNode >(const RCP< const CrsGraph<int, int, Kokkos::Compat::KokkosThreadsWrapperNode > > &graph);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraCrsGraphT<int, Kokkos::Compat::KokkosOpenMPWrapperNode >;
template RCP< const CrsGraph<int, int, Kokkos::Compat::KokkosOpenMPWrapperNode > > toXpetra<int, Kokkos::Compat::KokkosOpenMPWrapperNode>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<int, Kokkos::Compat::KokkosOpenMPWrapperNode >(const RCP< const CrsGraph<int, int, Kokkos::Compat::KokkosOpenMPWrapperNode > > &graph);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Kokkos::Compat::KokkosCudaWrapperNode default_node_type;
template class EpetraCrsGraphT<int, default_node_type >;
template RCP< const CrsGraph<int, int, default_node_type > > toXpetra<int, default_node_type>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<int, default_node_type >(const RCP< const CrsGraph<int, int, default_node_type > > &graph);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Kokkos::Compat::KokkosHIPWrapperNode default_node_type;
template class EpetraCrsGraphT<int, default_node_type >;
template RCP< const CrsGraph<int, int, default_node_type > > toXpetra<int, default_node_type>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<int, default_node_type >(const RCP< const CrsGraph<int, int, default_node_type > > &graph);
#endif
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
    (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
  template class EpetraCrsGraphT<long long, Xpetra::EpetraNode >;
  template RCP< const CrsGraph<int, long long, Xpetra::EpetraNode > > toXpetra<long long, Xpetra::EpetraNode>(const Epetra_CrsGraph &g);
  template const Epetra_CrsGraph & toEpetra<long long, Xpetra::EpetraNode >(const RCP< const CrsGraph<int, long long, Xpetra::EpetraNode > > &graph);
#endif
#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraCrsGraphT<long long, Kokkos::Compat::KokkosSerialWrapperNode >;
template RCP< const CrsGraph<int, long long, Kokkos::Compat::KokkosSerialWrapperNode > > toXpetra<long long, Kokkos::Compat::KokkosSerialWrapperNode>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<long long, Kokkos::Compat::KokkosSerialWrapperNode >(const RCP< const CrsGraph<int, long long, Kokkos::Compat::KokkosSerialWrapperNode > > &graph);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraCrsGraphT<long long, Kokkos::Compat::KokkosThreadsWrapperNode>;
template RCP< const CrsGraph<int, long long, Kokkos::Compat::KokkosThreadsWrapperNode > > toXpetra<long long, Kokkos::Compat::KokkosThreadsWrapperNode>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<long long, Kokkos::Compat::KokkosThreadsWrapperNode >(const RCP< const CrsGraph<int, long long, Kokkos::Compat::KokkosThreadsWrapperNode > > &graph);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraCrsGraphT<long long, Kokkos::Compat::KokkosOpenMPWrapperNode >;
template RCP< const CrsGraph<int, long long, Kokkos::Compat::KokkosOpenMPWrapperNode > > toXpetra<long long, Kokkos::Compat::KokkosOpenMPWrapperNode>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<long long, Kokkos::Compat::KokkosOpenMPWrapperNode >(const RCP< const CrsGraph<int, long long, Kokkos::Compat::KokkosOpenMPWrapperNode > > &graph);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Kokkos::Compat::KokkosCudaWrapperNode default_node_type;
template class EpetraCrsGraphT<long long, default_node_type >;
template RCP< const CrsGraph<int, long long, default_node_type > > toXpetra<long long, default_node_type>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<long long, default_node_type >(const RCP< const CrsGraph<int, long long, default_node_type > > &graph);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Kokkos::Compat::KokkosHIPWrapperNode default_node_type;
template class EpetraCrsGraphT<long long, default_node_type >;
template RCP< const CrsGraph<int, long long, default_node_type > > toXpetra<long long, default_node_type>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<long long, default_node_type >(const RCP< const CrsGraph<int, long long, default_node_type > > &graph);
#endif
#endif

} // namespace Xpetra
