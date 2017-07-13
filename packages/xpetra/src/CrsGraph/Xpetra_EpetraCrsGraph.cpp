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

#if 0
  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::EpetraCrsGraphT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra<EpetraGlobalOrdinal,Node>(rowMap), maxNumEntriesPerRow, toEpetra(pftype)))) { }
#endif
  // TODO: convert array size_t to int
  //   template<class EpetraGlobalOrdinal>
  //   EpetraCrsGraphT<EpetraGlobalOrdinal>::EpetraCrsGraphT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype)
  //     : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), NumEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)))) { }
#if 0
  template<class EpetraGlobalOrdinal, class Node>
  EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::EpetraCrsGraphT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra<EpetraGlobalOrdinal,Node>(rowMap), toEpetra<EpetraGlobalOrdinal,Node>(colMap), maxNumEntriesPerRow, toEpetra(pftype)))) { }
#endif
  // TODO: convert array size_t to int
  //   template<class EpetraGlobalOrdinal>
  //   EpetraCrsGraphT<EpetraGlobalOrdinal>::EpetraCrsGraphT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype)
  //     : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra<EpetraGlobalOrdinal,Node>(rowMap), toEpetra<EpetraGlobalOrdinal,Node>(colMap), NumEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)))) { }
#if 0
  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &indices) {
    XPETRA_MONITOR("EpetraCrsGraphT::insertGlobalIndices");

    GlobalOrdinal* indices_rawPtr = const_cast<GlobalOrdinal*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
    XPETRA_ERR_CHECK(graph_->InsertGlobalIndices(globalRow, indices.size(), indices_rawPtr));
  }
#endif

#if 0
  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::insertLocalIndices(int localRow, const ArrayView<const int> &indices) {
    XPETRA_MONITOR("EpetraCrsGraphT::insertLocalIndices");

    int* indices_rawPtr = const_cast<int*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
    XPETRA_ERR_CHECK(graph_->InsertMyIndices(localRow, indices.size(), indices_rawPtr));
  }
#endif
#if 0
  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &Indices) const {
    XPETRA_MONITOR("EpetraCrsGraphT::getGlobalRowView");

    int      numEntries;
    GlobalOrdinal    * eIndices;

    XPETRA_ERR_CHECK(graph_->ExtractGlobalRowView(GlobalRow, numEntries, eIndices));
    if (numEntries == 0) { eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    Indices = ArrayView<const GlobalOrdinal>(eIndices, numEntries);
  }
#endif
#if 0
  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::getLocalRowView(int LocalRow, ArrayView<const int> &indices) const {
    XPETRA_MONITOR("EpetraCrsGraphT::getLocalRowView");

    int      numEntries;
    int    * eIndices;

    XPETRA_ERR_CHECK(graph_->ExtractMyRowView(LocalRow, numEntries, eIndices));
    if (numEntries == 0) { eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
  }
#endif
#if 0
  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const RCP< ParameterList > &params){
    XPETRA_MONITOR("EpetraCrsGraphT::fillComplete");

    graph_->FillComplete(toEpetra<EpetraGlobalOrdinal,Node>(domainMap), toEpetra<EpetraGlobalOrdinal,Node>(rangeMap));
    bool doOptimizeStorage = true;
    if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
    if (doOptimizeStorage) graph_->OptimizeStorage();
  }
#endif
#if 0
  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::fillComplete(const RCP< ParameterList > &params) {
    XPETRA_MONITOR("EpetraCrsGraphT::fillComplete");

    graph_->FillComplete();
    bool doOptimizeStorage = true;
    if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
    if (doOptimizeStorage) graph_->OptimizeStorage();
  }
#endif
#if 0
  template<class EpetraGlobalOrdinal, class Node>
  std::string EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::description() const { XPETRA_MONITOR("EpetraCrsGraphT::description"); return "NotImplemented"; }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    XPETRA_MONITOR("EpetraCrsGraphT::describe");

    out << "EpetraCrsGraphT::describe : Warning, verbosity level is ignored by this method." << std::endl;
    const Epetra_BlockMap rowmap = graph_->RowMap();
    if (rowmap.Comm().MyPID() == 0) out << "** EpetraCrsGraphT **\n\nrowmap" << std::endl;
    rowmap.Print(out);
    graph_->Print(out);
  }
#endif

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
  //
#if 0
  // TODO: use toEpetra()
  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal, Node>::doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                                 const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraphT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraphT<GlobalOrdinal XPETRA_COMMA Node>, source, tSource, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraCrsGraphT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal XPETRA_COMMA Node>, importer, tImporter, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsGraph> v = tSource.getEpetra_CrsGraph();
    int err = graph_->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }
#endif
#if 0
  template<class EpetraGlobalOrdinal,class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal,Node>::doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                 const Import<LocalOrdinal, GlobalOrdinal, Node>& importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraphT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraphT<GlobalOrdinal XPETRA_COMMA Node>, dest, tDest, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraCrsGraphT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal XPETRA_COMMA Node>, importer, tImporter, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsGraph> v = tDest.getEpetra_CrsGraph();
    int err = graph_->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }
#endif
#if 0
  template<class EpetraGlobalOrdinal, class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal,Node>::doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                                 const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraphT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraphT<GlobalOrdinal XPETRA_COMMA Node>, source, tSource, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraCrsGraphT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal XPETRA_COMMA Node>, exporter, tExporter, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsGraph> v = tSource.getEpetra_CrsGraph();
    int err = graph_->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");

  }
#endif
#if 0
  template<class EpetraGlobalOrdinal,class Node>
  void EpetraCrsGraphT<EpetraGlobalOrdinal,Node>::doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                 const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraphT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraphT<GlobalOrdinal XPETRA_COMMA Node>, dest, tDest, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraCrsGraphT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal XPETRA_COMMA Node>, exporter, tExporter, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsGraph> v = tDest.getEpetra_CrsGraph();
    int err = graph_->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }
#endif

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
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
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraCrsGraphT<int, default_node_type >;
template RCP< const CrsGraph<int, int, default_node_type > > toXpetra<int, default_node_type>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<int, default_node_type >(const RCP< const CrsGraph<int, int, default_node_type > > &graph);
#endif // HAVE_XPETRA_TPETRA
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
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
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraCrsGraphT<long long, default_node_type >;
template RCP< const CrsGraph<int, long long, default_node_type > > toXpetra<long long, default_node_type>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<long long, default_node_type >(const RCP< const CrsGraph<int, long long, default_node_type > > &graph);
#endif // HAVE_XPETRA_TPETRA
#endif

} // namespace Xpetra
