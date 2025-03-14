// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_ISORROPIA_CRSGRAPH_H
#define TPETRA_ISORROPIA_CRSGRAPH_H

#include <EEP_Tpetra_Transform.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Tpetra {

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class CrsGraph;

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Map;


///
/** Generates an Tpetra::CrsGraph based on the repartitioning algorithms of Isorropia
 */
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Isorropia_CrsGraph : public StructuralSameTypeTransform< Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > {

  // Parameter list to specify partitioning algorithm for Isorropia.
  Teuchos::ParameterList partitionList_;

  // New graph for the repartitioned matrix. "String added by EEP, for testing"
  Teuchos::RCP< Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > NewGraph_;

 public:

  ///
  /* Destructor
   */
  ~Isorropia_CrsGraph();

  ///
  /* Constructor
   * input param part_method - type of Isorropia partitioning to use
   */
  Isorropia_CrsGraph( Teuchos::ParameterList& partitionList )
  : partitionList_(partitionList)
  {
    std::cout << "EEP Passing through trilinoscouplings/src/tpetra/EEP_Tpetra_Isorropia_CrsGraph.h Isorropia_CrsGraph::constructor()..." << std::endl;
  }

  ///
  /* Generates the Isorropia partitioned Tpetra::CrsGraph from the input object.
   */
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace Tpetra

#endif // TPETRA_ISORROPIA_CRSGRAPH_H
