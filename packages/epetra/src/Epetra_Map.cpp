
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
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
Epetra_Map::Epetra_Map(int NumGlobalElements, int NumMyElements, int * MyGlobalElements, 
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
Epetra_Map::~Epetra_Map(void) {}

//=============================================================================
Epetra_Map & Epetra_Map::operator= (const Epetra_Map& map) {
	if(this != &map)
		Epetra_BlockMap::operator=(map); // call this->Epetra_BlockMap::operator=
	return(*this);
}
