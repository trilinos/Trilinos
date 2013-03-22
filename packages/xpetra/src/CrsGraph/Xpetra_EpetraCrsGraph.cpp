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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
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
  const Epetra_CrsGraph & toEpetra(const RCP< const CrsGraph<int, int> > &graph) {
    XPETRA_RCP_DYNAMIC_CAST(const EpetraCrsGraph, graph, epetraGraph, "toEpetra");
    return *(epetraGraph->getEpetra_CrsGraph());
  }

  EpetraCrsGraph::EpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), maxNumEntriesPerRow, toEpetra(pftype)))) { }

  // TODO: convert array size_t to int
  //   EpetraCrsGraph::EpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype)
  //     : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), NumEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)))) { }

  EpetraCrsGraph::EpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), toEpetra(colMap), maxNumEntriesPerRow, toEpetra(pftype)))) { }

  // TODO: convert array size_t to int
  //   EpetraCrsGraph::EpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype)
  //     : graph_(Teuchos::rcp(new Epetra_CrsGraph(Copy, toEpetra(rowMap), toEpetra(colMap), NumEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)))) { }

  void EpetraCrsGraph::insertGlobalIndices(int globalRow, const ArrayView<const int> &indices) {
    XPETRA_MONITOR("EpetraCrsGraph::insertGlobalIndices");

    int* indices_rawPtr = const_cast<int*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
    XPETRA_ERR_CHECK(graph_->InsertGlobalIndices(globalRow, indices.size(), indices_rawPtr));
  }

  void EpetraCrsGraph::insertLocalIndices(int localRow, const ArrayView<const int> &indices) {
    XPETRA_MONITOR("EpetraCrsGraph::insertLocalIndices");

    int* indices_rawPtr = const_cast<int*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
    XPETRA_ERR_CHECK(graph_->InsertMyIndices(localRow, indices.size(), indices_rawPtr));
  }

  void EpetraCrsGraph::getGlobalRowView(int GlobalRow, ArrayView<const int> &Indices) const {
    XPETRA_MONITOR("EpetraCrsGraph::getGlobalRowView");

    int      numEntries;
    int    * eIndices;

    XPETRA_ERR_CHECK(graph_->ExtractGlobalRowView(GlobalRow, numEntries, eIndices));
    if (numEntries == 0) { eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    Indices = ArrayView<const int>(eIndices, numEntries);
  }

  void EpetraCrsGraph::getLocalRowView(int LocalRow, ArrayView<const int> &indices) const {
    XPETRA_MONITOR("EpetraCrsGraph::getLocalRowView");

    int      numEntries;
    int    * eIndices;

    XPETRA_ERR_CHECK(graph_->ExtractMyRowView(LocalRow, numEntries, eIndices));
    if (numEntries == 0) { eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
  }

  void EpetraCrsGraph::fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const RCP< ParameterList > &params){
    XPETRA_MONITOR("EpetraCrsGraph::fillComplete");

    graph_->FillComplete(toEpetra(domainMap), toEpetra(rangeMap));
    bool doOptimizeStorage = true;
    if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
    if (doOptimizeStorage) graph_->OptimizeStorage();
  }

  void EpetraCrsGraph::fillComplete(const RCP< ParameterList > &params) {
    XPETRA_MONITOR("EpetraCrsGraph::fillComplete");

    graph_->FillComplete();
    bool doOptimizeStorage = true;
    if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
    if (doOptimizeStorage) graph_->OptimizeStorage();
  }

  std::string EpetraCrsGraph::description() const { XPETRA_MONITOR("EpetraCrsGraph::description"); return "NotImplemented"; }

  void EpetraCrsGraph::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    XPETRA_MONITOR("EpetraCrsGraph::describe");

    out << "EpetraCrsGraph::describe : Warning, verbosity level is ignored by this method." << std::endl;
    const Epetra_BlockMap rowmap = graph_->RowMap();
    if (rowmap.Comm().MyPID() == 0) out << "** EpetraCrsGraph **\n\nrowmap" << std::endl;
    rowmap.Print(out);
    graph_->Print(out);
  }

  // TODO: move that elsewhere
  RCP< const CrsGraph<int, int> > toXpetra(const Epetra_CrsGraph &g) {
    RCP<const Epetra_CrsGraph> const_graph = rcp(new Epetra_CrsGraph(g));

    RCP<Epetra_CrsGraph> graph = Teuchos::rcp_const_cast<Epetra_CrsGraph>(const_graph); //TODO: can I avoid the const_cast ?
    return rcp( new Xpetra::EpetraCrsGraph(graph) );
  }
  //

  // TODO: use toEpetra()
  void EpetraCrsGraph::doImport(const DistObject<int, int, int> &source,
                                 const Import<int, int> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraph::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraph, source, tSource, "Xpetra::EpetraCrsGraph::doImport only accept Xpetra::EpetraCrsGraph as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraCrsGraph::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsGraph> v = tSource.getEpetra_CrsGraph();
    int err = graph_->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraCrsGraph::doExport(const DistObject<int, int, int> &dest,
                                 const Import<int, int>& importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraph::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraph, dest, tDest, "Xpetra::EpetraCrsGraph::doImport only accept Xpetra::EpetraCrsGraph as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraCrsGraph::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsGraph> v = tDest.getEpetra_CrsGraph();
    int err = graph_->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraCrsGraph::doImport(const DistObject<int, int, int> &source,
                                 const Export<int, int>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraph::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraph, source, tSource, "Xpetra::EpetraCrsGraph::doImport only accept Xpetra::EpetraCrsGraph as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraCrsGraph::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsGraph> v = tSource.getEpetra_CrsGraph();
    int err = graph_->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");

  }

  void EpetraCrsGraph::doExport(const DistObject<int, int, int> &dest,
                                 const Export<int, int>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsGraph::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsGraph, dest, tDest, "Xpetra::EpetraCrsGraph::doImport only accept Xpetra::EpetraCrsGraph as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraCrsGraph::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsGraph> v = tDest.getEpetra_CrsGraph();
    int err = graph_->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

} // namespace Xpetra
