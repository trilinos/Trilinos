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

  using tCrsGraph = Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;

  // Parameter list to specify partitioning algorithm for Isorropia.
  Teuchos::ParameterList partitionList_;

  // New graph for the repartitioned matrix.
  Teuchos::RCP< tCrsGraph > NewGraph_;

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
  typename Isorropia_CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::NewTypeRef operator()( typename Isorropia_CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::OriginalTypeRef orig );

};

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Isorropia_CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
~Isorropia_CrsGraph()
{
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename Isorropia_CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::NewTypeRef
Isorropia_CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
operator()( typename Isorropia_CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::OriginalTypeRef orig )
{
  std::cout << "EEP Entering trilinoscouplings/src/tpetra/EEP_Tpetra_Isorropia_CrsGraph.cpp Isorropia_CrsGraph::operator()..." << std::endl;

  this->origObj_ = &orig;

#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES // AquiToDo
#else
  std::cout << "EEP In trilinoscouplings/src/tpetra/EEP_Tpetra_Isorropia_CrsGraph.cpp Isorropia_CrsGraph::operator(), pos 001" << std::endl;
  if  (sizeof(GlobalOrdinal) == 8) {
    // There is nothing we can do since Isorropia does not support 64-bit integers, just create a copy of the original graph.
    NewGraph_ = Teuchos::rcp( new tCrsGraph( orig ) );
  } 
  else
#endif
  if (orig.getGlobalNumRows() == 0) {
    std::cout << "EEP In trilinoscouplings/src/tpetra/EEP_Tpetra_Isorropia_CrsGraph.cpp Isorropia_CrsGraph::operator(), pos 002" << std::endl;
    // If there is nothing to do, just create a copy of the original empty graph.
    NewGraph_ = Teuchos::rcp( new tCrsGraph( orig ) );
  }
  else {
    std::cout << "EEP In trilinoscouplings/src/tpetra/EEP_Tpetra_Isorropia_CrsGraph.cpp Isorropia_CrsGraph::operator(), pos 003" << std::endl; // Aqui
    try {
      //NewGraph_ = Teuchos::rcp( Isorropia::Tpetra::createBalancedCopy( orig, partitionList_) ); // AquiToDo
    }
    catch(std::exception& e) {
      std::cout << "Isorropia::Tpetra::createBalancedCopy threw exception '" << e.what() << "' on proc " << orig.getComm()->getRank() << std::endl;
    }
  }

  std::cout << "EEP In trilinoscouplings/src/tpetra/EEP_Tpetra_Isorropia_CrsGraph.cpp Isorropia_CrsGraph::operator(), pos 004" << std::endl;
  // Set the raw pointer to the new graph.
  // This should be OK, since the destructor does not attempt to destroy the raw pointer.
  this->newObj_ = NewGraph_.get();

  std::cout << "EEP Leaving trilinoscouplings/src/tpetra/EEP_Tpetra_Isorropia_CrsGraph.cpp Isorropia_CrsGraph::operator()" << std::endl;
  return *NewGraph_;
}

} //namespace Tpetra

#endif // TPETRA_ISORROPIA_CRSGRAPH_H
