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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_XPETRA_SUP_GRAPH_UTILS_HPP_
#define PACKAGES_XPETRA_SUP_GRAPH_UTILS_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_MatrixMatrix.hpp"

namespace Xpetra {

/*!
  @class GraphUtils
  @brief Xpetra utility class for common graph-related routines

  The routines should be independent from Epetra/Tpetra and be purely implemented in Xpetra.

*/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class GraphUtils {
#undef XPETRA_GRAPHUTILS_SHORT
#include "Xpetra_UseShortNames.hpp"

public:

  /*!
    \@brief Get the exponential of a CrsGraph

    @param graph Graph to calculate the exponential of
    @param levelFill Parameter to control the recursive calculation of the exponential

    @return Exponential CrsGraph with increased fill-in
  */
  static RCP<const CrsGraph> getExponentialCrsGraph(const RCP<const CrsGraph>& graph, int levelFill, Teuchos::FancyOStream &fos) {

    RCP<Matrix> adjacency = MatrixFactory::Build(graph);
    adjacency->setAllToScalar(1.0);
    adjacency->fillComplete();

    RCP<Matrix> exponentialAdjacency = adjacency;

    for(int power=levelFill; power>1; power--)
      exponentialAdjacency = MatrixMatrix::Multiply(*exponentialAdjacency, false, *adjacency, false, fos);

    return exponentialAdjacency->getCrsGraph();
  }

};

} // end namespace Xpetra

#define XPETRA_GRAPHUTILS_SHORT

#endif // PACKAGES_XPETRA_SUP_GRAPH_UTILS_HPP_
