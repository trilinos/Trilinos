
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_Map.h"


//==============================================================================
// Epetra_Map constructor for a Epetra-defined uniform linear distribution of elements.
Epetra_Map::Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_BlockMap( NumGlobalElements, 1, IndexBase, Comm) // Map is just a special case of BlockMap
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
Epetra_Map::~Epetra_Map(void) 
{
}
