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

#ifndef GALERI_CARTESIAN2D_H
#define GALERI_CARTESIAN2D_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

namespace Galeri {

namespace Maps {

inline
Epetra_Map* 
Cartesian2D(const Epetra_Comm& Comm, const int nx, const int ny,
            const int mx, const int my)
{
  if (nx <= 0 || ny <= 0 || mx <= 0 || my <= 0 || (mx > nx) || (my > ny))
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Cartesian2D()",
                    "nx = " + toString(nx) +
                    ", ny = " + toString(ny) +
                    ", mx = " + toString(mx) +
                    ", my = " + toString(my)));

  // how to divide the axis
  int modx = (nx + (nx % mx)) / mx;
  int mody = (ny + (ny % my)) / my;

  int MyPID = Comm.MyPID(), startx, starty, endx, endy;
  int xpid = MyPID % mx;
  int ypid = MyPID / mx;

  startx = xpid * modx; 
  if ((xpid + 1) * modx < nx) endx = (xpid + 1) * modx;
  else                        endx = nx;
  starty = ypid * mody; 
  if ((ypid + 1) * mody < ny) endy = (ypid + 1) * mody;
  else                        endy = ny;

  int NumMyElements = (endx - startx) * (endy - starty);
  vector<int> MyGlobalElements(NumMyElements);
  int count = 0;

  for (int i = startx ; i < endx ; ++i) 
    for (int j = starty ; j < endy ; ++j) 
      MyGlobalElements[count++] = i + j * nx;

  return(new Epetra_Map(nx * ny,  NumMyElements, &MyGlobalElements[0],
                        0, Comm));

} // Cartesian2D()

} // namespace BlockMaps
} // namespace Galeri
#endif
