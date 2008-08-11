
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
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
// ************************************************************************
//@HEADER

#include "Epetra_Map.h"


//==============================================================================
// Epetra_Map constructor for a Epetra-defined uniform linear distribution of elements.
Epetra_Map::Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_BlockMap(NumGlobalElements, 1, IndexBase, Comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
//==============================================================================
// Epetra_Map constructor for a user-defined linear distribution of constant block size elements.
Epetra_Map::Epetra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_BlockMap(NumGlobalElements, NumMyElements, 1, IndexBase, Comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
//==============================================================================
// Epetra_Map constructor for a user-defined arbitrary distribution of constant block size elements.
Epetra_Map::Epetra_Map(int NumGlobalElements, int NumMyElements,
                       const int * MyGlobalElements,
                       int IndexBase, const Epetra_Comm& Comm)
  : Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, 1, IndexBase, Comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
//==============================================================================
Epetra_Map::Epetra_Map(const Epetra_Map& map)
  : Epetra_BlockMap(map) // Map is just a special case of BlockMap
{
}

//==============================================================================
Epetra_Map::~Epetra_Map(void)
{
}

//=============================================================================
Epetra_Map & Epetra_Map::operator= (const Epetra_Map& map)
{
  if(this != &map) {
    Epetra_BlockMap::operator=(map); // call this->Epetra_BlockMap::operator=
  }
  return(*this);
}
