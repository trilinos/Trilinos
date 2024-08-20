// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_RANDOM_H
#define GALERI_RANDOM_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Util.h"
#include <limits>

namespace Galeri {
namespace Maps {

template<typename int_type>
inline
Epetra_Map* 
TRandom(const Epetra_Comm& Comm, const int_type n)
{
  if (n <= 0 || n > std::numeric_limits<int>::max())
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Random()",
                    "n = " + toString(n)));
  
  // This is the idea: I create the map on proc 0, then I broadcast
  // it to all procs. This is not very efficient, but saves some MPI calls.
      
  std::vector<int> part(n);
      
  if (Comm.MyPID() == 0) 
  {
    Epetra_Util Util;

    for (int_type i = 0 ; i < n ; ++i) 
    {
      unsigned int r = Util.RandomInt();	
      part[i] = r%(Comm.NumProc());
    }
  }
      
  Comm.Broadcast(&part[0], n, 0);
      
  // count the elements assigned to this proc
  int NumMyElements = 0;
  for (int_type i = 0 ; i < n; ++i) 
    if (part[i] == Comm.MyPID()) NumMyElements++;

  // get the loc2global list
  std::vector<int_type> MyGlobalElements(NumMyElements);
  int count = 0;
  for (int_type i = 0 ; i < n ; ++i)
    if (part[i] == Comm.MyPID()) 
      MyGlobalElements[count++] = i;

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  return(new Epetra_Map(n, NumMyElements, (NumMyElements ? &MyGlobalElements[0] : NULL), 0, Comm));
#else
  return(new Epetra_Map);
#endif

} // TRandom()

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map* 
Random(const Epetra_Comm& Comm, const int n) {
  return TRandom<int>(Comm, n);
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map* 
Random64(const Epetra_Comm& Comm, const long long n) {
  return TRandom<long long>(Comm, n);
}
#endif

} // namespace BlockMaps
} // namespace Galeri
#endif
