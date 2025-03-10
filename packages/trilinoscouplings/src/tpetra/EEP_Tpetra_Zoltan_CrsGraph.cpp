// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <EEP_Tpetra_Zoltan_CrsGraph.h>
#include <Isorropia_Epetra.hpp>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

namespace EpetraExt {

Zoltan_CrsGraph::
~Zoltan_CrsGraph()
{
}

Zoltan_CrsGraph::NewTypeRef
Zoltan_CrsGraph::
operator()( OriginalTypeRef orig )
{
  std::cout << "EEP Entering trilinoscouplings/src/tpetra/EEP_Tpetra_Zoltan_CrsGraph.cpp Zoltan_CrsGraph::operator()..." << std::endl;
  origObj_ = &orig;

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  std::cout << "EEP In trilinoscouplings/src/tpetra/EEP_Tpetra_Zoltan_CrsGraph.cpp Zoltan_CrsGraph::operator(), pos 001" << std::endl;
  const Epetra_BlockMap & OldMap = orig.RowMap();
  if(OldMap.GlobalIndicesLongLong()) {
    // There is nothing we can do since Isorropia does not support Epetra 64-bit integers, just create a copy of the original graph.
    NewGraph_ = Teuchos::rcp( new Epetra_CrsGraph( orig ) );
  } 
  else
#endif
  std::cout << "EEP In trilinoscouplings/src/tpetra/EEP_Tpetra_Zoltan_CrsGraph.cpp Zoltan_CrsGraph::operator(), pos 002" << std::endl;
  if (orig.NumGlobalRows() == 0)
  {
    std::cout << "EEP In trilinoscouplings/src/tpetra/EEP_Tpetra_Zoltan_CrsGraph.cpp Zoltan_CrsGraph::operator(), pos 003" << std::endl;
    // If there is nothing to do, just create a copy of the original empty graph.
    NewGraph_ = Teuchos::rcp( new Epetra_CrsGraph( orig ) );
  }
  else {
    std::cout << "EEP In trilinoscouplings/src/tpetra/EEP_Tpetra_Zoltan_CrsGraph.cpp Zoltan_CrsGraph::operator(), pos 004" << std::endl; // Aqui
    try {
      NewGraph_ =
        Teuchos::rcp( Isorropia::Epetra::createBalancedCopy( orig, partitionList_) );
    }
    catch(std::exception& e) {
      std::cout << "Isorropia::create_balanced_copy threw exception '" << e.what() << "' on proc " << orig.Comm().MyPID() << std::endl;
    }
  }

  std::cout << "EEP In trilinoscouplings/src/tpetra/EEP_Tpetra_Zoltan_CrsGraph.cpp Zoltan_CrsGraph::operator(), pos 005" << std::endl;
  // Set the raw pointer to the new graph.
  // This should be OK, since the destructor does not attempt to destroy the raw pointer.
  newObj_ = NewGraph_.get();

  std::cout << "EEP Leaving trilinoscouplings/src/tpetra/EEP_Tpetra_Zoltan_CrsGraph.cpp Zoltan_CrsGraph::operator()" << std::endl;
  return *NewGraph_;
}

} // namespace Tpetra

