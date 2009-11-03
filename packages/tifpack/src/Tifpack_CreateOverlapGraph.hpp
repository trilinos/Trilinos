/*@HEADER
// ***********************************************************************
// 
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef TIFPACK_CREATEOVERLAPGRAPH_HPP
#define TIFPACK_CREATEOVERLAPGRAPH_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace Tifpack {

//! Constructs an overlapped graph for use with Tifpack preconditioners.
/**
  If OverlapLevel is 0, then the overlapped graph is the input_graph.
 */
template<class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > CreateOverlapGraph(const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& input_graph, int OverlapLevel)
{
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> GraphType;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> MapType;
  typedef Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> ImportType;

  TEST_FOR_EXCEPTION(OverlapLevel < 0, std::runtime_error, "Tifpack::CreateOverlapGraph: OverlapLevel must be >= 0.");

  Teuchos::RCP<GraphType> OverlapGraph = input_graph;

  const int numProcs = input_graph->getMap()->getComm()->getSize();
  if (OverlapLevel == 0 || numProcs < 2) return OverlapGraph;

  Teuchos::RCP<const MapType> OverlapRowMap = input_graph->getRowMap();

  Teuchos::RCP<GraphType> OldGraph;
  Teuchos::RCP<const MapType> OldRowMap;
  const Teuchos::RCP<const MapType> DomainMap = input_graph->getDomainMap();
  const Teuchos::RCP<const MapType> RangeMap = input_graph->getRangeMap();

  for (int level=1; level <= OverlapLevel; level++) {
    OldGraph = OverlapGraph;
    OldRowMap = OverlapRowMap;

    Teuchos::RCP<const ImportType> OverlapImporter = OldGraph->getImporter();
    OverlapRowMap = OverlapImporter->getTargetMap();
    if (level<OverlapLevel) {
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

}//namespace Tifpack

#endif // TIFPACK_CREATEOVERLAPGRAPH_HPP
