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
  template<class GlobalOrdinal>
  const Epetra_CrsGraph & toEpetra(const RCP< const CrsGraph<int, GlobalOrdinal> > &graph) {
    XPETRA_RCP_DYNAMIC_CAST(const EpetraCrsGraphT<GlobalOrdinal>, graph, epetraGraph, "toEpetra");
    return *(epetraGraph->getEpetra_CrsGraph());
  }

  template<class EpetraGlobalOrdinal>
  EpetraCrsGraphT<EpetraGlobalOrdinal>::EpetraCrsGraphT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), maxNumEntriesPerRow, toEpetra(pftype)))) { }

  // TODO: convert array size_t to int
  //   template<class EpetraGlobalOrdinal>
  //   EpetraCrsGraphT<EpetraGlobalOrdinal>::EpetraCrsGraphT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype)
  //     : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), NumEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)))) { }

  template<class EpetraGlobalOrdinal>
  EpetraCrsGraphT<EpetraGlobalOrdinal>::EpetraCrsGraphT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), toEpetra(colMap), maxNumEntriesPerRow, toEpetra(pftype)))) { }

  // TODO: convert array size_t to int
  //   template<class EpetraGlobalOrdinal>
  //   EpetraCrsGraphT<EpetraGlobalOrdinal>::EpetraCrsGraphT(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype)
  //     : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), toEpetra(colMap), NumEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)))) { }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &indices) {
    XPETRA_MONITOR("EpetraCrsGraphT::insertGlobalIndices");

    GlobalOrdinal* indices_rawPtr = const_cast<GlobalOrdinal*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
    XPETRA_ERR_CHECK(graph_->InsertGlobalIndices(globalRow, indices.size(), indices_rawPtr));
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::insertLocalIndices(int localRow, const ArrayView<const int> &indices) {
    XPETRA_MONITOR("EpetraCrsGraphT::insertLocalIndices");

    int* indices_rawPtr = const_cast<int*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
    XPETRA_ERR_CHECK(graph_->InsertMyIndices(localRow, indices.size(), indices_rawPtr));
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &Indices) const {
    XPETRA_MONITOR("EpetraCrsGraphT::getGlobalRowView");

    int      numEntries;
    GlobalOrdinal    * eIndices;

    XPETRA_ERR_CHECK(graph_->ExtractGlobalRowView(GlobalRow, numEntries, eIndices));
    if (numEntries == 0) { eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    Indices = ArrayView<const GlobalOrdinal>(eIndices, numEntries);
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::getLocalRowView(int LocalRow, ArrayView<const int> &indices) const {
    XPETRA_MONITOR("EpetraCrsGraphT::getLocalRowView");

    int      numEntries;
    int    * eIndices;

    XPETRA_ERR_CHECK(graph_->ExtractMyRowView(LocalRow, numEntries, eIndices));
    if (numEntries == 0) { eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const RCP< ParameterList > &params){
    XPETRA_MONITOR("EpetraCrsGraphT::fillComplete");

    graph_->FillComplete(toEpetra(domainMap), toEpetra(rangeMap));
    bool doOptimizeStorage = true;
    if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
    if (doOptimizeStorage) graph_->OptimizeStorage();
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::fillComplete(const RCP< ParameterList > &params) {
    XPETRA_MONITOR("EpetraCrsGraphT::fillComplete");

    graph_->FillComplete();
    bool doOptimizeStorage = true;
    if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
    if (doOptimizeStorage) graph_->OptimizeStorage();
  }

  template<class EpetraGlobalOrdinal>
  std::string EpetraCrsGraphT<EpetraGlobalOrdinal>::description() const { XPETRA_MONITOR("EpetraCrsGraphT::description"); return "NotImplemented"; }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    XPETRA_MONITOR("EpetraCrsGraphT::describe");

    out << "EpetraCrsGraphT::describe : Warning, verbosity level is ignored by this method." << std::endl;
    const Epetra_BlockMap rowmap = graph_->RowMap();
    if (rowmap.Comm().MyPID() == 0) out << "** EpetraCrsGraphT **\n\nrowmap" << std::endl;
    rowmap.Print(out);
    graph_->Print(out);
  }

  // TODO: move that elsewhere
  template<class GlobalOrdinal>
  RCP<const CrsGraph<int, GlobalOrdinal> >
  toXpetra (const Epetra_CrsGraph &g)
  {
    RCP<const Epetra_CrsGraph> const_graph = rcp (new Epetra_CrsGraph (g));
    RCP<Epetra_CrsGraph> graph =
      Teuchos::rcp_const_cast<Epetra_CrsGraph> (const_graph);
    return rcp (new Xpetra::EpetraCrsGraphT<GlobalOrdinal>(graph));
  }
  //

  // TODO: use toEpetra()
  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::doImport(const DistObject<GlobalOrdinal, int, GlobalOrdinal> &source,
                                 const Import<int, GlobalOrdinal> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraphT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraphT<GlobalOrdinal>, source, tSource, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraCrsGraphT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, importer, tImporter, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsGraph> v = tSource.getEpetra_CrsGraph();
    int err = graph_->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::doExport(const DistObject<GlobalOrdinal, int, GlobalOrdinal> &dest,
                                 const Import<int, GlobalOrdinal>& importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraphT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraphT<GlobalOrdinal>, dest, tDest, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraCrsGraphT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, importer, tImporter, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsGraph> v = tDest.getEpetra_CrsGraph();
    int err = graph_->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::doImport(const DistObject<GlobalOrdinal, int, GlobalOrdinal> &source,
                                 const Export<int, GlobalOrdinal>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraphT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraphT<GlobalOrdinal>, source, tSource, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraCrsGraphT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal>, exporter, tExporter, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsGraph> v = tSource.getEpetra_CrsGraph();
    int err = graph_->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");

  }

  template<class EpetraGlobalOrdinal>
  void EpetraCrsGraphT<EpetraGlobalOrdinal>::doExport(const DistObject<GlobalOrdinal, int, GlobalOrdinal> &dest,
                                 const Export<int, GlobalOrdinal>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraphT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraphT<GlobalOrdinal>, dest, tDest, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraCrsGraphT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal>, exporter, tExporter, "Xpetra::EpetraCrsGraphT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<const Epetra_CrsGraph> v = tDest.getEpetra_CrsGraph();
    int err = graph_->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template class EpetraCrsGraphT<int>;
template RCP< const CrsGraph<int, int> > toXpetra<int>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<int>(const RCP< const CrsGraph<int, int> > &graph);
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template class EpetraCrsGraphT<long long>;
template RCP< const CrsGraph<int, long long> > toXpetra<long long>(const Epetra_CrsGraph &g);
template const Epetra_CrsGraph & toEpetra<long long>(const RCP< const CrsGraph<int, long long> > &graph);
#endif

} // namespace Xpetra
