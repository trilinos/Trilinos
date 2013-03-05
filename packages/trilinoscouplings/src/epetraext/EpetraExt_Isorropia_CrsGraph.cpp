// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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

