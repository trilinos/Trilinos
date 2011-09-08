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

#ifndef GALERI_NODECARTESIAN2D_H
#define GALERI_NODECARTESIAN2D_H
#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#ifndef HAVE_MPI
#include "Epetra_Comm.h"
#else
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"

namespace Galeri {

namespace Maps {

inline
Epetra_Map* 
NodeCartesian2D(const Epetra_Comm& Comm, const Epetra_Comm & NodeComm, const int MyNodeID,
                const int nx, const int ny,
                const int ndx, const int ndy, const int px, const int py )
{
  if (nx <= 0 || ny <= 0 || (ndx > nx) || (ndy > ny) || px <= 0 || py <= 0)
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::NodeCartesian2D()",
                    "nx = " + toString(nx) +
                    ", ny = " + toString(ny) +
                    ", ndx = " + toString(ndx) +
                    ", ndy = " + toString(ndy) +
                    ", px = " + toString(px) +
                    ", py = " + toString(py)));
  int NodePID=NodeComm.MyPID();
  // Compute nodal position
  int nxpid = MyNodeID % ndx;
  int nypid = MyNodeID / ndx;

  // Compute position w/i node
  int pxpid = NodePID % px;
  int pypid = NodePID / px;

  // Compute effective PID, mx, my
  int mx = ndx*px;
  int my = ndy*py;
  int ppn=px*py;
  int EffectivePID = nypid*ndx*ppn + nxpid*px + pypid*mx + pxpid;

  // Build the map based on the Effective PID - code copied from Galeri_NodeCartesian2D.h
  int xpid = EffectivePID % mx;
  int ypid = EffectivePID / mx;

  int PerProcSmallXDir = (int) (((double) nx)/((double) mx));
  int PerProcSmallYDir = (int) (((double) ny)/((double) my));
  int NBigXDir         = nx - PerProcSmallXDir*mx;
  int NBigYDir         = ny - PerProcSmallYDir*my;
  int startx,starty,endx,endy;
  
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
  vector<int> MyGlobalElements(NumMyElements);
  int count = 0;

  for (int i = startx ; i < endx ; ++i) 
    for (int j = starty ; j < endy ; ++j) 
      MyGlobalElements[count++] = i + j * nx;

  return(new Epetra_Map(nx * ny,  NumMyElements, &MyGlobalElements[0],
                        0, Comm));

  
} // NodeCartesian2D()

} // namespace BlockMaps
} // namespace Galeri
#endif
