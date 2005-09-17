#ifndef GALERI_INTERLACED_H
#define GALERI_INTERLACED_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

namespace Galeri {
namespace Maps {

inline
Epetra_Map* 
Interlaced(Epetra_Comm& Comm, int NumGlobalElements)
{
  if (NumGlobalElements <= 0)
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Interlaced()",
                    "n = " + toString(NumGlobalElements)));
                    
  // this is a funky map. Nodes are assigned so that
  // node 0 is given to proc 0, node 1 to proc 1, and
  // node i to proc i%NumProcs. Probably not the best, but it
  // results in decompositions with lots of boundary nodes.

  int NumProcs = Comm.NumProc();
  int MyPID    = Comm.MyPID();

  int NumMyElements = NumGlobalElements / NumProcs;
  if (MyPID < NumGlobalElements % NumProcs) NumMyElements++;

  int count = 0;
  vector<int> MyGlobalElements(NumMyElements);

  for (int i = 0 ; i < NumGlobalElements ; ++i) 
  {
    if (i%NumProcs == MyPID) 
      MyGlobalElements[count++] = i;
  }

  if (count != NumMyElements)
    throw(Exception(__FILE__, __LINE__,
                    "Something went wrong in Maps::Interlaced()",
                    "count = " + toString(count) 
                    + ", NumMyElements = " + toString(NumMyElements)));

  return(new Epetra_Map(NumGlobalElements, NumMyElements,
                          &MyGlobalElements[0], 0, Comm));
} // Interlaced()

} // namespace Maps
} // namespace Galeri
#endif
