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
//  Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#ifndef GALERI_CARTESIAN3D_H
#define GALERI_CARTESIAN3D_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

namespace Galeri {
namespace Maps {

inline
Epetra_Map* 
Cartesian3D(const Epetra_Comm& Comm, const int nx, const int ny, const int nz,
            const int mx, const int my, const int mz)
{
  if (nx <= 0 || ny <= 0 || nz <= 0 ||
      mx <= 0 || my <= 0 || mz <= 0 ||
      (mx > nx) || (my > ny) || (mz > nz))
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Cartesian3D()",
                    "nx = " + toString(nx) +
                    ", ny = " + toString(ny) +
                    ", nz = " + toString(nz) +
                    ", mx = " + toString(mx) +
                    ", my = " + toString(my) +
                    ", mz = " + toString(mz)));

  int MyPID = Comm.MyPID();
  int startx, starty, startz, endx, endy, endz;

  int modx = (nx +(nx %mx )) / mx;
  int mody = (ny +(ny %my )) / my;
  int modz = (nz +(nz %mz )) / mz;

  int mxy  = mx * my;
  int zpid = MyPID / mxy;
  int xpid = (MyPID % mxy) % mx;
  int ypid = (MyPID % mxy) / mx;

  startx = xpid * modx;
  if ((xpid + 1) * modx < nx) endx = (xpid + 1) * modx;
  else                        endx = nx;
  starty = ypid * mody;
  if ((ypid + 1) * mody < ny) endy = (ypid + 1) * mody;
  else                        endy = ny;
  startz = zpid * modz;
  if ((zpid + 1) * modz < nz) endz = (zpid + 1) * modz;
  else                        endz = nz;

  int NumMyElements = (endx - startx) * (endy - starty) * (endz - startz);
  vector<int> MyGlobalElements(NumMyElements);
  int count = 0;

  for (int i = startx ; i < endx ; ++i) 
    for (int j = starty ; j < endy ; ++j) 
      for (int k = startz ; k < endz ; ++k) 
        MyGlobalElements[count++] = i + j * nx +k * (nx * ny);

  return(new Epetra_Map (-1, NumMyElements, &MyGlobalElements[0], 0, Comm));

} // Cartesian3D()

} // namespace Maps
} // namespace Galeri
#endif
