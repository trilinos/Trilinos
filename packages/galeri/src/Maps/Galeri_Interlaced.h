// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
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
