/*@HEADER
// ***********************************************************************
// 
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
template<class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > CreateOverlapGraph(const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& input_graph, int OverlapLevel)
{
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> GraphType;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> MapType;
  typedef Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> ImportType;

  TEUCHOS_TEST_FOR_EXCEPTION(OverlapLevel < 0, std::runtime_error, "Ifpack2::CreateOverlapGraph: OverlapLevel must be >= 0.");

  Teuchos::RCP<GraphType> OverlapGraph;

  const int numProcs = input_graph->getMap()->getComm()->getSize();
  if (OverlapLevel == 0 || numProcs < 2) return input_graph;

  Teuchos::RCP<const MapType> OverlapRowMap = input_graph->getRowMap();

  Teuchos::RCP<const GraphType> OldGraph;
  Teuchos::RCP<const MapType> OldRowMap;
  const Teuchos::RCP<const MapType> DomainMap = input_graph->getDomainMap();
  const Teuchos::RCP<const MapType> RangeMap = input_graph->getRangeMap();

  for (int level=0; level < OverlapLevel; level++) {
    OldGraph = OverlapGraph;
    OldRowMap = OverlapRowMap;

    Teuchos::RCP<const ImportType> OverlapImporter; 
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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > CreateOverlapMatrix(const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& input_graph, int OverlapLevel)
{
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> MatrixType;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> MapType;
  typedef Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> ImportType;

  TEUCHOS_TEST_FOR_EXCEPTION(OverlapLevel < 0, std::runtime_error, "Ifpack2::CreateOverlapMatrix: OverlapLevel must be >= 0.");

  Teuchos::RCP<MatrixType> OverlapGraph;

  const int numProcs = input_graph->getMap()->getComm()->getSize();
  if (OverlapLevel == 0 || numProcs < 2) return input_graph;

  Teuchos::RCP<const MapType> OverlapRowMap = input_graph->getRowMap();

  Teuchos::RCP<const MatrixType> OldGraph;
  Teuchos::RCP<const MapType> OldRowMap;
  const Teuchos::RCP<const MapType> DomainMap = input_graph->getDomainMap();
  const Teuchos::RCP<const MapType> RangeMap = input_graph->getRangeMap();

  for (int level=0; level < OverlapLevel; level++) {
    OldGraph = OverlapGraph;
    OldRowMap = OverlapRowMap;

    Teuchos::RCP<const ImportType> OverlapImporter; 
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
