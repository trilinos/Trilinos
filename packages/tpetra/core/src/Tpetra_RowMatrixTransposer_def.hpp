// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_ROWMATRIXTRANSPOSER_DEF_HPP
#define TPETRA_ROWMATRIXTRANSPOSER_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_shortSort.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_SortCrs.hpp"

namespace Tpetra {

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
RowMatrixTransposer (const Teuchos::RCP<const crs_matrix_type>& origMatrix,
                     const std::string& label)
  : origMatrix_ (origMatrix), label_ (label)
{}

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
createTranspose (const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  using Teuchos::RCP;
  // Do the local transpose
  RCP<crs_matrix_type> transMatrixWithSharedRows = createTransposeLocal (params);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  const std::string prefix = std::string ("Tpetra ") + label_ + ": ";
  using Teuchos::TimeMonitor;
  TimeMonitor MM (*TimeMonitor::getNewTimer (prefix + "Transpose TAFC"));
#endif

  // If transMatrixWithSharedRows has an exporter, that's what we
  // want.  If it doesn't, the rows aren't actually shared, and we're
  // done!
  using export_type = Export<LocalOrdinal, GlobalOrdinal, Node>;
  RCP<const export_type> exporter =
    transMatrixWithSharedRows->getGraph ()->getExporter ();
  if (exporter.is_null ()) {
    return transMatrixWithSharedRows;
  }
  else {
    Teuchos::ParameterList labelList;
#ifdef HAVE_TPETRA_MMM_TIMINGS
    labelList.set("Timer Label", label_);
#endif
    if(! params.is_null ()) {
      const char paramName[] = "compute global constants";
      labelList.set (paramName, params->get (paramName, true));
    }
    // Use the Export object to do a fused Export and fillComplete.
    // This always sorts the local matrix after communication, so
    //   no need to set "sorted = false" in parameters.
    return exportAndFillCompleteCrsMatrix<crs_matrix_type>
      (transMatrixWithSharedRows, *exporter, Teuchos::null,
       Teuchos::null, Teuchos::rcpFromRef (labelList));
  }
}


// mfh 03 Feb 2013: In a definition outside the class like this, the
// return value is considered outside the class scope (for things like
// resolving typedefs), but the arguments are considered inside the
// class scope.
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
createTransposeLocal (const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using import_type = Tpetra::Import<LO, GO, Node>;
  using export_type = Tpetra::Export<LO, GO, Node>;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix = std::string("Tpetra ") + label_ + ": ";
  using Teuchos::TimeMonitor;
  TimeMonitor MM (*TimeMonitor::getNewTimer (prefix + "Transpose Local"));
#endif

  const bool sort = [&] () {
    constexpr bool sortDefault = true; // see #4607 discussion
    const char sortParamName[] = "sort";
    return params.get () == nullptr ? sortDefault :
      params->get (sortParamName, sortDefault);
  } ();

  const LO lclNumRows (origMatrix_->getLocalNumRows ());

  RCP<const crs_matrix_type> crsMatrix =
    rcp_dynamic_cast<const crs_matrix_type> (origMatrix_);
  if (crsMatrix.is_null ()) {
    auto rowMap = origMatrix_->getRowMap ();
    if (rowMap->isOneToOne ()) {
      Teuchos::Array<size_t> numEntPerRow (lclNumRows);
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        numEntPerRow[lclRow] = origMatrix_->getNumEntriesInLocalRow (lclRow);
      }
      auto colMap = origMatrix_->getColMap ();

      RCP<crs_matrix_type> crsMatrix_nc =
        rcp (new crs_matrix_type (rowMap, colMap, numEntPerRow ()));

      // When source & target Maps are same, Import just copies.
      import_type imp (rowMap, rowMap);
      crsMatrix_nc->doImport (*origMatrix_, imp, Tpetra::REPLACE);
      crsMatrix_nc->fillComplete (origMatrix_->getDomainMap (),
                                  origMatrix_->getRangeMap ());
      crsMatrix = crsMatrix_nc;
    }
    else {
      TEUCHOS_ASSERT( false ); // not implemented (it wasn't before)
    }
  }

  using local_matrix_device_type = typename crs_matrix_type::local_matrix_device_type;

  local_matrix_device_type lclMatrix = crsMatrix->getLocalMatrixDevice ();
  local_matrix_device_type lclTransposeMatrix = KokkosSparse::Impl::transpose_matrix(lclMatrix);
  if (sort)
    KokkosSparse::sort_crs_matrix(lclTransposeMatrix);

  // Prebuild the importers and exporters the no-communication way,
  // flipping the importers and exporters around.
  const auto origExport = origMatrix_->getGraph ()->getExporter ();
  RCP<const import_type> myImport = origExport.is_null () ?
    Teuchos::null : rcp (new import_type (*origExport));
  const auto origImport = origMatrix_->getGraph ()->getImporter ();
  RCP<const export_type> myExport = origImport.is_null () ?
    Teuchos::null : rcp (new export_type (*origImport));

  RCP<Teuchos::ParameterList> graphParams = Teuchos::null;
  if(!sort) {
    graphParams = rcp(new Teuchos::ParameterList);
    graphParams->set("sorted", false);
  }

  return rcp (new crs_matrix_type (lclTransposeMatrix,
                                   origMatrix_->getColMap (),
                                   origMatrix_->getRowMap (),
                                   origMatrix_->getRangeMap (),
                                   origMatrix_->getDomainMap (),
                                   myImport, myExport, graphParams));
}

/*************************************************************************/

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
BlockCrsMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
BlockCrsMatrixTransposer (const Teuchos::RCP<const bcrs_matrix_type>& origMatrix,
                     const std::string& label)
  : origMatrix_ (origMatrix), label_ (label)
{}

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
Teuchos::RCP<BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
BlockCrsMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
createTranspose (const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  using Teuchos::RCP;
  // Do the local transpose
  RCP<bcrs_matrix_type> transMatrixWithSharedRows = createTransposeLocal (params);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  const std::string prefix = std::string ("Tpetra ") + label_ + ": ";
  using Teuchos::TimeMonitor;
  TimeMonitor MM (*TimeMonitor::getNewTimer (prefix + "Transpose TAFC"));
#endif

  // If transMatrixWithSharedRows has an exporter, that's what we
  // want.  If it doesn't, the rows aren't actually shared, and we're
  // done!
  using export_type = Export<LocalOrdinal, GlobalOrdinal, Node>;
  RCP<const export_type> exporter =
    transMatrixWithSharedRows->getGraph ()->getExporter ();
  if (exporter.is_null ()) {
    return transMatrixWithSharedRows;
  }
  else {
    Teuchos::ParameterList labelList;
#ifdef HAVE_TPETRA_MMM_TIMINGS
    labelList.set("Timer Label", label_);
#endif
    if(! params.is_null ()) {
      const char paramName[] = "compute global constants";
      labelList.set (paramName, params->get (paramName, true));
    }
    // Use the Export object to do a fused Export and fillComplete.
    // This always sorts the local matrix after communication, so
    //   no need to set "sorted = false" in parameters.
    return exportAndFillCompleteBlockCrsMatrix<bcrs_matrix_type>
      (transMatrixWithSharedRows, *exporter);
  }
}


// mfh 03 Feb 2013: In a definition outside the class like this, the
// return value is considered outside the class scope (for things like
// resolving typedefs), but the arguments are considered inside the
// class scope.
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
Teuchos::RCP<BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
BlockCrsMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
createTransposeLocal (const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using import_type = Tpetra::Import<LO, GO, Node>;
  using export_type = Tpetra::Export<LO, GO, Node>;
  using crs_graph_type = typename bcrs_matrix_type::crs_graph_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix = std::string("Tpetra ") + label_ + ": ";
  using Teuchos::TimeMonitor;
  TimeMonitor MM (*TimeMonitor::getNewTimer (prefix + "Transpose Local"));
#endif

  RCP<const bcrs_matrix_type> crsMatrix =
    rcp_dynamic_cast<const bcrs_matrix_type> (origMatrix_);

  if(crsMatrix.is_null()) 
    TEUCHOS_ASSERT( false ); // not implemented 

  using local_matrix_device_type = typename bcrs_matrix_type::local_matrix_device_type;

  typename local_matrix_device_type::values_type values ;
  RCP<const crs_graph_type> graph;
  {
    local_matrix_device_type lclMatrix = crsMatrix->getLocalMatrixDevice ();

    local_matrix_device_type lclTransposeMatrix = KokkosSparse::Impl::transpose_bsr_matrix(lclMatrix);    

    // BlockCrs requires that we sort stuff
    KokkosSparse::sort_crs_matrix(lclTransposeMatrix);
    values = lclTransposeMatrix.values;
    
    // Prebuild the importers and exporters the no-communication way,
    // flipping the importers and exporters around.
    const auto origExport = origMatrix_->getGraph ()->getExporter ();
    RCP<const import_type> myImport = origExport.is_null () ?
      Teuchos::null : rcp (new import_type (*origExport));
    const auto origImport = origMatrix_->getGraph ()->getImporter ();
    RCP<const export_type> myExport = origImport.is_null () ?
      Teuchos::null : rcp (new export_type (*origImport));

    RCP<Teuchos::ParameterList> graphParams = Teuchos::null;
    
    // Make the Transpose Graph
    graph = rcp(new crs_graph_type(lclTransposeMatrix.graph,                                   
                                   origMatrix_->getColMap (),
                                   origMatrix_->getRowMap (),
                                   origMatrix_->getGraph()->getRangeMap (),
                                   origMatrix_->getGraph()->getDomainMap (),
                                   myImport,
                                   myExport,
                                   graphParams));
  }
  // Now make the matrix
  return rcp (new bcrs_matrix_type (*graph,
                                    values,
                                    origMatrix_->getBlockSize()));
}
//


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIXTRANSPOSER_INSTANT(SCALAR,LO,GO,NODE) \
  template class RowMatrixTransposer< SCALAR, LO , GO , NODE >;\
  template class BlockCrsMatrixTransposer< SCALAR, LO , GO , NODE >;

} // namespace Tpetra

#endif
