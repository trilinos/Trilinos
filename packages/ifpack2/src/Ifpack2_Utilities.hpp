// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_UTILITIES_HPP
#define IFPACK2_UTILITIES_HPP

#include "Ifpack2_ConfigDefs.hpp"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsGraph.hpp"

//! @file Ifpack2_Utilities.hpp
//! @brief File for utility functions.

namespace Ifpack2 {

namespace Details {


  /*! @brief Compute and return the graph of the diagonal of the input graph.

      This is copied from Tpetra::Experimental::BlockCrsMatrix. 
  */
  template <class graph_type>
  Teuchos::RCP<Tpetra::CrsGraph<typename graph_type::local_ordinal_type, typename graph_type::global_ordinal_type, typename graph_type::node_type> >
  computeDiagonalGraph (graph_type const &graph)
  {
    typedef typename graph_type::local_ordinal_type  LO;
    typedef typename graph_type::global_ordinal_type GO;
    typedef typename graph_type::node_type           NO;
    typedef Tpetra::Map<LO, GO, NO> map_type;
    typedef Tpetra::CrsGraph<LO, GO, NO> crs_graph_type;

    const size_t maxDiagEntPerRow = 1;
    // NOTE (mfh 12 Aug 2014) We could also pass in the column Map
    // here.  However, we still would have to do LID->GID lookups to
    // make sure that we are using the correct diagonal column
    // indices, so it probably wouldn't help much.
    Teuchos::RCP<graph_type> diagonalGraph;
    diagonalGraph = Teuchos::rcp(new crs_graph_type(graph.getRowMap(), maxDiagEntPerRow));
    const map_type& meshRowMap = *(graph.getRowMap());

    Teuchos::Array<GO> diagGblColInds(maxDiagEntPerRow);

    for (LO lclRowInd = meshRowMap.getMinLocalIndex(); lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd) {
      const GO gblRowInd = meshRowMap.getGlobalElement(lclRowInd);
      diagGblColInds[0] = gblRowInd;
      diagonalGraph->insertGlobalIndices(gblRowInd, diagGblColInds());
    }

    diagonalGraph->fillComplete(graph.getDomainMap(), graph.getRangeMap());
    return diagonalGraph;
  }


  //! Transform to canonical form of preconditioner name.
  std::string canonicalize(const std::string& precType);


}// namespace Details

}// namespace Ifpack2

#endif //IFPACK2_UTILITIES_HPP
