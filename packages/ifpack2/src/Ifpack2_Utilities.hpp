/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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
    diagonalGraph = Teuchos::rcp(new crs_graph_type(graph.getRowMap(), maxDiagEntPerRow, Tpetra::StaticProfile));
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
