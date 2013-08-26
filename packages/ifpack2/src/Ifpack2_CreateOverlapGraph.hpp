/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_CREATEOVERLAPGRAPH_HPP
#define IFPACK2_CREATEOVERLAPGRAPH_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace Ifpack2 {

//! Constructs an overlapped graph for use with Ifpack2 preconditioners.
/**
  If OverlapLevel is 0, then the overlapped graph is the input_graph.
 */
template<class GraphType>
Teuchos::RCP<const GraphType> CreateOverlapGraph(const Teuchos::RCP<const GraphType>& input_graph, int OverlapLevel)
{
  typedef typename GraphType::map_type map_type;
  typedef Tpetra::Import<typename GraphType::local_ordinal_type,typename GraphType::global_ordinal_type,typename GraphType::node_type> import_type;

  TEUCHOS_TEST_FOR_EXCEPTION(OverlapLevel < 0, std::runtime_error, "Ifpack2::CreateOverlapGraph: OverlapLevel must be >= 0.");

  Teuchos::RCP<GraphType> OverlapGraph;

  const int numProcs = input_graph->getMap()->getComm()->getSize();
  if (OverlapLevel == 0 || numProcs < 2) return input_graph;

  Teuchos::RCP<const map_type> OverlapRowMap = input_graph->getRowMap();

  Teuchos::RCP<const GraphType> OldGraph;
  Teuchos::RCP<const map_type> OldRowMap;
  const Teuchos::RCP<const map_type> DomainMap = input_graph->getDomainMap();
  const Teuchos::RCP<const map_type> RangeMap = input_graph->getRangeMap();

  for (int level=0; level < OverlapLevel; level++) {
    OldGraph = OverlapGraph;
    OldRowMap = OverlapRowMap;

    Teuchos::RCP<const import_type> OverlapImporter; 
    if(level==0) OverlapImporter = input_graph->getImporter();
    else OverlapImporter = OldGraph->getImporter();

    OverlapRowMap = OverlapImporter->getTargetMap();
    if (level<OverlapLevel-1) {
      OverlapGraph = Teuchos::rcp( new GraphType(OverlapRowMap, 0) );
    }
    else {
      // On last iteration, we want to filter out all columns except those that
      // correspond to rows in the graph.  This assures that our matrix is square
      OverlapGraph = Teuchos::rcp( new GraphType(OverlapRowMap, OverlapRowMap, 0) );
    }

    OverlapGraph->doImport(*input_graph, *OverlapImporter, Tpetra::INSERT);
    OverlapGraph->fillComplete(DomainMap, RangeMap);
  }

  return OverlapGraph;
}

template<class MatrixType>
Teuchos::RCP<const MatrixType> CreateOverlapMatrix(const Teuchos::RCP<const MatrixType>& input_graph, int OverlapLevel)
{
  typedef typename MatrixType::map_type map_type;
  typedef Tpetra::Import<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> import_type;

  TEUCHOS_TEST_FOR_EXCEPTION(OverlapLevel < 0, std::runtime_error, "Ifpack2::CreateOverlapMatrix: OverlapLevel must be >= 0.");

  Teuchos::RCP<MatrixType> OverlapGraph;

  const int numProcs = input_graph->getMap()->getComm()->getSize();
  if (OverlapLevel == 0 || numProcs < 2) return input_graph;

  Teuchos::RCP<const map_type> OverlapRowMap = input_graph->getRowMap();

  Teuchos::RCP<const MatrixType> OldGraph;
  Teuchos::RCP<const map_type> OldRowMap;
  const Teuchos::RCP<const map_type> DomainMap = input_graph->getDomainMap();
  const Teuchos::RCP<const map_type> RangeMap = input_graph->getRangeMap();

  for (int level=0; level < OverlapLevel; level++) {
    OldGraph = OverlapGraph;
    OldRowMap = OverlapRowMap;

    Teuchos::RCP<const import_type> OverlapImporter; 
    if(level==0) OverlapImporter = input_graph->getGraph()->getImporter();
    else OverlapImporter = OldGraph->getGraph()->getImporter();

    OverlapRowMap = OverlapImporter->getTargetMap();
    if (level<OverlapLevel-1) {
      OverlapGraph = Teuchos::rcp( new MatrixType(OverlapRowMap, 0) );
    }
    else {
      // On last iteration, we want to filter out all columns except those that
      // correspond to rows in the graph.  This assures that our matrix is square
      OverlapGraph = Teuchos::rcp( new MatrixType(OverlapRowMap, OverlapRowMap, 0) );
    }

    OverlapGraph->doImport(*input_graph, *OverlapImporter, Tpetra::INSERT);
    OverlapGraph->fillComplete(DomainMap, RangeMap);
  }

  return OverlapGraph;
}



}//namespace Ifpack2

#endif // IFPACK2_CREATEOVERLAPGRAPH_HPP
