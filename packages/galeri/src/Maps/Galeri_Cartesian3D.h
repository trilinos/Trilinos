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

  int mxy  = mx * my;
  int zpid = MyPID / mxy;
  int xpid = (MyPID % mxy) % mx;
  int ypid = (MyPID % mxy) / mx;

  int PerProcSmallXDir = (int) (((double) nx)/((double) mx));
  int PerProcSmallYDir = (int) (((double) ny)/((double) my));
  int PerProcSmallZDir = (int) (((double) nz)/((double) mz));
  int NBigXDir         = nx - PerProcSmallXDir*mx;
  int NBigYDir         = ny - PerProcSmallYDir*my;
  int NBigZDir         = nz - PerProcSmallZDir*mz;

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

  if (zpid < NBigZDir) startz =  zpid*(PerProcSmallZDir+1);
  else                 startz = (zpid-NBigZDir)*PerProcSmallZDir +
                                      NBigZDir*(PerProcSmallZDir+1);
  endz = startz + PerProcSmallZDir;
  if ( zpid < NBigZDir) endz++;

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
