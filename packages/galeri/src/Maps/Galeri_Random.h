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

#ifndef GALERI_RANDOM_H
#define GALERI_RANDOM_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Util.h"

namespace Galeri {
namespace Maps {

inline
Epetra_Map* 
Random(const Epetra_Comm& Comm, const int n)
{
  if (n <= 0)
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Random()",
                    "n = " + toString(n)));
  
  // This is the idea: I create the map on proc 0, then I broadcast
  // it to all procs. This is not very efficient, but saves some MPI calls.
      
  vector<int> part(n);
      
  if (Comm.MyPID() == 0) 
  {
    Epetra_Util Util;

    for (int i = 0 ; i < n ; ++i) 
    {
      unsigned int r = Util.RandomInt();	
      part[i] = r%(Comm.NumProc());
    }
  }
      
  Comm.Broadcast(&part[0], n, 0);
      
  // count the elements assigned to this proc
  int NumMyElements = 0;
  for (int i = 0 ; i < n; ++i) 
    if (part[i] == Comm.MyPID()) NumMyElements++;

  // get the loc2global list
  vector<int> MyGlobalElements(NumMyElements);
  int count = 0;
  for (int i = 0 ; i < n ; ++i)
    if (part[i] == Comm.MyPID()) 
      MyGlobalElements[count++] = i;

  return(new Epetra_Map(n, NumMyElements, &MyGlobalElements[0], 0, Comm));

} // Random()

} // namespace BlockMaps
} // namespace Galeri
#endif
