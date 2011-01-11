// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Piro_Epetra_StokhosNOXObserver.hpp"

Piro::Epetra::StokhosNOXObserver::StokhosNOXObserver (
     Teuchos::RCP<NOX::Epetra::Observer> noxObserver_,
     const Epetra_Map& map_,
     const int sz_) :
  noxObserver(noxObserver_),
  map(map_),
  numSGBlocks(sz_)
{
 //if (noxObserver == Teuchos::null) cout << "XXX1" << endl;
}

void Piro::Epetra::StokhosNOXObserver::observeSolution(
    const Epetra_Vector& solution)
{

  // Copy into block vector, so Block access is available
  EpetraExt::BlockVector blockVec(View, map, solution);

  Teuchos::RCP<Epetra_Vector> deterministicSizedSolnVec;

  for (int i=0; i< numSGBlocks; i++) {
    deterministicSizedSolnVec = blockVec.GetBlock(i);

    if (noxObserver != Teuchos::null)
      noxObserver->observeSolution(*deterministicSizedSolnVec);
  }

}
