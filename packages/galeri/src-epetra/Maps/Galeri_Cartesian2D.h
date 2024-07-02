// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_CARTESIAN2D_H
#define GALERI_CARTESIAN2D_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

namespace Galeri {

namespace Maps {

template<typename int_type>
inline
Epetra_Map* 
TCartesian2D(const Epetra_Comm& Comm, const int nx, const int ny,
            const int mx, const int my)
{
  if (nx <= 0 || ny <= 0 || mx <= 0 || my <= 0 || (mx > nx) || (my > ny))
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Cartesian2D()",
                    "nx = " + toString(nx) +
                    ", ny = " + toString(ny) +
                    ", mx = " + toString(mx) +
                    ", my = " + toString(my)));

  int MyPID = Comm.MyPID(), startx, starty, endx, endy;
  int xpid = MyPID % mx;
  int ypid = MyPID / mx;

  int PerProcSmallXDir = (int) (((double) nx)/((double) mx));
  int PerProcSmallYDir = (int) (((double) ny)/((double) my));
  int NBigXDir         = nx - PerProcSmallXDir*mx;
  int NBigYDir         = ny - PerProcSmallYDir*my;

  if (xpid < NBigXDir) startx =  xpid*(PerProcSmallXDir+1);
  else                 startx = (xpid-NBigXDir)*PerProcSmallXDir + 
                                      NBigXDir*(PerProcSmallXDir+1);
  endx = startx + PerProcSmallXDir;
  if (xpid < NBigXDir) endx++;

  if (ypid < NBigYDir) starty =  ypid*(PerProcSmallYDir+1);
  else                 starty = (ypid-NBigYDir)*PerProcSmallYDir +
                                      NBigYDir*(PerProcSmallYDir+1);
  endy = starty + PerProcSmallYDir;
  if ( ypid < NBigYDir) endy++;

  int NumMyElements = (endx - startx) * (endy - starty);
  std::vector<int_type> MyGlobalElements(NumMyElements);
  int count = 0;

  for (int i = startx ; i < endx ; ++i) 
    for (int j = starty ; j < endy ; ++j) 
      MyGlobalElements[count++] = i + j * nx;

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  return(new Epetra_Map(nx * ny,  NumMyElements, &MyGlobalElements[0],
                        0, Comm));
#else
  return(new Epetra_Map);
#endif

} // TCartesian2D()

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map* 
Cartesian2D(const Epetra_Comm& Comm, const int nx, const int ny,
            const int mx, const int my) {
  return TCartesian2D<int>(Comm, nx, ny, mx, my);
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map* 
Cartesian2D64(const Epetra_Comm& Comm, const int nx, const int ny,
            const int mx, const int my) {
  return TCartesian2D<long long>(Comm, nx, ny, mx, my);
}
#endif

} // namespace Maps
} // namespace Galeri
#endif
