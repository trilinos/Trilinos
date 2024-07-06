// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_INTERLACED_H
#define GALERI_INTERLACED_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include <limits>

namespace Galeri {
namespace Maps {

template<typename int_type>
inline
Epetra_Map* 
TInterlaced(Epetra_Comm& Comm, int_type NumGlobalElements)
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

  if (NumGlobalElements / NumProcs > std::numeric_limits<int>::max())
    throw(Exception(__FILE__, __LINE__,
                    "Too high NumGlobalElements to Maps::Interlaced()",
                    "n = " + toString(NumGlobalElements)));

  int NumMyElements = (int) (NumGlobalElements / NumProcs);
  if (MyPID < NumGlobalElements % NumProcs) NumMyElements++;

  int count = 0;
  std::vector<int_type> MyGlobalElements(NumMyElements);

  for (int_type i = 0 ; i < NumGlobalElements ; ++i) 
  {
    if (i%NumProcs == MyPID) 
      MyGlobalElements[count++] = i;
  }

  if (count != NumMyElements)
    throw(Exception(__FILE__, __LINE__,
                    "Something went wrong in Maps::Interlaced()",
                    "count = " + toString(count) 
                    + ", NumMyElements = " + toString(NumMyElements)));

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  return(new Epetra_Map(NumGlobalElements, NumMyElements,
                          &MyGlobalElements[0], 0, Comm));
#else
  return(new Epetra_Map);
#endif
} // TInterlaced()

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map* 
Interlaced(Epetra_Comm& Comm, int NumGlobalElements) {
  return TInterlaced<int>(Comm, NumGlobalElements);
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map* 
Interlaced64(Epetra_Comm& Comm, long long NumGlobalElements) {
  return TInterlaced<long long>(Comm, NumGlobalElements);
}
#endif

} // namespace Maps
} // namespace Galeri
#endif
