// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <EpetraExt_Isorropia_CrsGraph.h>
#include <Isorropia_Epetra.hpp>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

namespace EpetraExt {

Isorropia_CrsGraph::
~Isorropia_CrsGraph()
{
}

Isorropia_CrsGraph::NewTypeRef
Isorropia_CrsGraph::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  const Epetra_BlockMap & OldMap = orig.RowMap();
  if(OldMap.GlobalIndicesLongLong()) {
    // There is nothing we can do since Isorropia does not support Epetra 64-bit integers, just create a copy of the original graph.
    NewGraph_ = Teuchos::rcp( new Epetra_CrsGraph( orig ) );
  } 
  else
#endif
  if (orig.NumGlobalRows() == 0)
  {
    // If there is nothing to do, just create a copy of the original empty graph.
    NewGraph_ = Teuchos::rcp( new Epetra_CrsGraph( orig ) );
  }
  else {
    try {
      NewGraph_ =
        Teuchos::rcp( Isorropia::Epetra::createBalancedCopy( orig, partitionList_) );
    }
    catch(std::exception& e) {
      std::cout << "Isorropia::create_balanced_copy threw exception '" << e.what() << "' on proc " << orig.Comm().MyPID() << std::endl;
    }
  }

  // Set the raw pointer to the new graph.
  // This should be OK, since the destructor does not attempt to destroy the raw pointer.
  newObj_ = NewGraph_.get();

  return *NewGraph_;
}

} // namespace EpetraExt

