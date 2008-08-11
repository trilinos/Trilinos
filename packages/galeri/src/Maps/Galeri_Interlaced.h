// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
