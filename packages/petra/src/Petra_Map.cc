
#include "Petra_Map.h"


//==============================================================================
// Petra_Map constructor for a Petra-defined uniform linear distribution of elements.
Petra_Map::Petra_Map(int NumGlobalElements, int IndexBase, const Petra_Comm& Comm)
  : Petra_BlockMap( NumGlobalElements, 1, IndexBase, Comm)
{

}
//==============================================================================
// Petra_Map constructor for a user-defined linear distribution of constant block size elements.
Petra_Map::Petra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Petra_Comm& Comm)
  : Petra_BlockMap(NumGlobalElements, NumMyElements, 1, IndexBase, Comm)
{
}
//==============================================================================
// Petra_Map constructor for a user-defined arbitrary distribution of constant block size elements.
Petra_Map::Petra_Map(int NumGlobalElements, int NumMyElements, int * MyGlobalElements, 
		     int IndexBase, const Petra_Comm& Comm)
  : Petra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, 1, IndexBase, Comm)
{
}
//==============================================================================
Petra_Map::Petra_Map(const Petra_Map& map)
  : Petra_BlockMap(map)
{
}

//==============================================================================
Petra_Map::~Petra_Map(void) 
{
}
